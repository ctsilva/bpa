#include <IO.h>
#include <algorithm>
#include <bpa.h>
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <glm/glm.hpp>
#include <iostream>
#include <map>
#include <numbers>
#include <set>

using namespace bpa;

namespace {
	// Points on the unit sphere, `n` of them in a Fibonacci spiral, with outward normals.
	auto fibonacciSphere(int n, glm::dvec3 center = {}, double radius = 1) -> std::vector<Point> {
		std::vector<Point> points;
		const auto golden = std::numbers::pi * (3 - std::sqrt(5.0));
		for (auto i = 0; i < n; i++) {
			const auto z = 1 - 2.0 * (i + 0.5) / n;
			const auto r = std::sqrt(1 - z * z);
			const auto phi = golden * i;
			const glm::dvec3 v{r * std::cos(phi), r * std::sin(phi), z};
			points.push_back({center + radius * v, v});
		}
		return points;
	}

	auto measuredReconstruct(const std::vector<Point>& points, double radius, const Options& options = {}) -> std::vector<Face> {
		const auto start = std::chrono::steady_clock::now();
		auto result = reconstruct(points, radius, options);
		const auto seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
		std::cerr << "Points: " << points.size() << ", Triangles: " << result.size() << ", T/s: " << result.size() / seconds << '\n';
		return result;
	}

	struct MeshStats {
		std::size_t boundaryEdges = 0; // edges with one triangle
		std::size_t badEdges = 0;      // edges with more than two triangles, or two with the same direction
		std::size_t components = 0;
		std::size_t verticesUsed = 0;
		long eulerCharacteristic = 0;
	};

	auto stats(const std::vector<Face>& faces, std::size_t nPoints) -> MeshStats {
		MeshStats s;
		std::map<std::pair<std::uint32_t, std::uint32_t>, int> directed; // half-edge -> count
		std::set<std::uint32_t> used;
		std::vector<std::uint32_t> parent(nPoints);
		for (std::uint32_t i = 0; i < nPoints; i++)
			parent[i] = i;
		const auto root = [&](std::uint32_t x) {
			while (parent[x] != x)
				x = parent[x] = parent[parent[x]];
			return x;
		};
		for (const auto& t : faces) {
			for (auto k = 0; k < 3; k++) {
				directed[{t[k], t[(k + 1) % 3]}]++;
				used.insert(t[k]);
			}
			parent[root(t[0])] = parent[root(t[1])] = root(t[2]);
		}
		std::set<std::pair<std::uint32_t, std::uint32_t>> undirected;
		for (const auto& [e, n] : directed) {
			undirected.insert({std::min(e.first, e.second), std::max(e.first, e.second)});
			const auto reverse = directed.count({e.second, e.first}) ? directed.at({e.second, e.first}) : 0;
			if (n > 1 || reverse > 1)
				s.badEdges++;
			else if (reverse == 0)
				s.boundaryEdges++;
		}
		std::set<std::uint32_t> roots;
		for (const auto& t : faces)
			roots.insert(root(t[0]));
		s.components = roots.size();
		s.verticesUsed = used.size();
		s.eulerCharacteristic = long(used.size()) - long(undirected.size()) + long(faces.size());
		return s;
	}
} // namespace

TEST_CASE("closed sphere", "[reconstruct]") {
	// 1.5 x the mean spacing: every point used, closed, Euler characteristic 2
	const auto cloud = fibonacciSphere(2000);
	const auto mesh = measuredReconstruct(cloud, 0.119);
	const auto s = stats(mesh, cloud.size());
	CHECK(mesh.size() == 2 * cloud.size() - 4);
	CHECK(s.verticesUsed == cloud.size());
	CHECK(s.boundaryEdges == 0);
	CHECK(s.badEdges == 0);
	CHECK(s.components == 1);
	CHECK(s.eulerCharacteristic == 2);
}

TEST_CASE("triangles face outward", "[reconstruct]") {
	const auto cloud = fibonacciSphere(500);
	const auto mesh = reconstruct(cloud, 0.24);
	REQUIRE(!mesh.empty());
	for (const auto& t : mesh) {
		const auto n = glm::cross(cloud[t[1]].pos - cloud[t[0]].pos, cloud[t[2]].pos - cloud[t[0]].pos);
		CHECK(glm::dot(n, cloud[t[0]].normal) > 0);
	}
}

TEST_CASE("two spheres need two seeds", "[reconstruct]") {
	auto cloud = fibonacciSphere(1000, {0, 0, 0});
	const auto second = fibonacciSphere(1000, {5, 0, 0});
	cloud.insert(cloud.end(), second.begin(), second.end());
	const auto mesh = measuredReconstruct(cloud, 0.17);
	const auto s = stats(mesh, cloud.size());
	CHECK(mesh.size() == 2 * (2 * 1000 - 4));
	CHECK(s.components == 2);
	CHECK(s.boundaryEdges == 0);
	CHECK(s.verticesUsed == cloud.size());
}

TEST_CASE("minComponent drops the small sphere", "[reconstruct]") {
	auto cloud = fibonacciSphere(1000, {0, 0, 0});
	const auto small = fibonacciSphere(50, {5, 0, 0}, 0.2);
	cloud.insert(cloud.end(), small.begin(), small.end());
	const auto all = reconstruct(cloud, 0.17);
	CHECK(stats(all, cloud.size()).components == 2);
	const auto large = reconstruct(cloud, 0.17, Options{.minComponent = 200});
	CHECK(stats(large, cloud.size()).components == 1);
	CHECK(large.size() == 2 * 1000 - 4);
}

namespace {
	// A jittered grid of nx x ny points in the plane z = 0 with normal +z, spacing h.
	auto planePatch(int nx, int ny, double h, glm::dvec3 origin, unsigned seed) -> std::vector<Point> {
		std::vector<Point> points;
		auto state = seed;
		const auto jitter = [&] {
			state = state * 1664525u + 1013904223u;
			return (state / 4294967296.0 - 0.5) * 0.4 * h;
		};
		for (auto j = 0; j < ny; j++)
			for (auto i = 0; i < nx; i++)
				points.push_back({origin + glm::dvec3{i * h + jitter(), j * h + jitter(), 0}, {0, 0, 1}});
		return points;
	}
} // namespace

TEST_CASE("multiple radii on uneven sampling", "[reconstruct]") {
	// a dense patch next to a patch sampled 2.5 times as coarsely
	const auto h = 1.0 / 40;
	auto cloud = planePatch(20, 40, h, {0, 0, 0}, 3);
	const auto sparse = planePatch(8, 16, 2.5 * h, {20 * h, 0, 0}, 4);
	cloud.insert(cloud.end(), sparse.begin(), sparse.end());

	// one small radius: the sparse patch is not reached
	const auto one = reconstruct(cloud, 1.3 * h);
	const auto s1 = stats(one, cloud.size());
	CHECK(s1.verticesUsed < cloud.size());
	CHECK(s1.badEdges == 0);

	// two radii: everything is covered by one disk; the first pass is the single-radius mesh
	const auto two = reconstruct(cloud, std::vector<double>{1.3 * h, 1.3 * 2.5 * h});
	const auto s2 = stats(two, cloud.size());
	CHECK(s2.verticesUsed == cloud.size());
	CHECK(s2.components == 1);
	CHECK(s2.badEdges == 0);
	CHECK(s2.eulerCharacteristic == 1);
	CHECK(two.size() > one.size());
	CHECK(std::equal(one.begin(), one.end(), two.begin()));

	// the order of the radii does not matter
	CHECK(reconstruct(cloud, std::vector<double>{1.3 * 2.5 * h, 1.3 * h}) == two);
}

TEST_CASE("radius too small for the spacing", "[reconstruct]") {
	const auto cloud = fibonacciSphere(2000);
	CHECK(reconstruct(cloud, 0.01).empty());
	CHECK_THROWS(reconstruct(cloud, std::vector<double>{}));
	CHECK_THROWS(reconstruct(cloud, -1.0));
}

TEST_CASE("bunny", "[reconstruct]") {
	const auto cloud = loadXYZ(TEST_DATA_DIR "/bunny.xyz");
	const auto mesh = measuredReconstruct(cloud, 0.002);
	const auto s = stats(mesh, cloud.size());
	CHECK(mesh.size() > 60000);
	CHECK(s.badEdges == 0);
	saveSTL(TEST_OUTPUT_DIR "/bunny_mesh.stl", cloud, mesh); // to look at, next to the test binary
}
