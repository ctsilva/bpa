#include <IO.h>
#include <bpa.h>
#include <catch2/catch_test_macros.hpp>
#include <chrono>
#include <glm/glm.hpp>
#include <iostream>
#include <numbers>

using namespace bpa;

namespace {
	auto createSphericalCloud(int slices, int stacks) -> std::vector<Point> {
		std::vector<Point> points;
		points.emplace_back(Point{{0, 0, -1}, {0, 0, -1}});
		for (auto slice = 0; slice < slices; slice++) {
			for (auto stack = 1; stack < stacks; stack++) {
				const auto yaw = (static_cast<double>(slice) / slices) * 2 * std::numbers::pi;
				const auto z = std::sin((static_cast<double>(stack) / stacks - 0.5) * std::numbers::pi);
				const auto r = std::sqrt(1 - z * z);

				glm::dvec3 v;
				v.x = static_cast<double>(r * std::sin(yaw));
				v.y = static_cast<double>(r * std::cos(yaw));
				v.z = static_cast<double>(z);
				points.push_back({v, glm::normalize(v - glm::dvec3{})});
			}
		}
		points.emplace_back(Point{{0, 0, 1}, {0, 0, 1}});
		return points;
	}

	auto measuredReconstruct(const std::vector<Point>& points, double radius) -> std::vector<Face> {
		const auto start = std::chrono::high_resolution_clock::now();
		auto result = reconstruct(points, radius);
		const auto end = std::chrono::high_resolution_clock::now();
		const auto seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
		std::cerr << "Points: " << points.size() << ", Triangles: " << result.size() << ", T/s: " << result.size() / seconds << '\n';
		return result;
	}
} // namespace

TEST_CASE("sphere_36_18", "[reconstruct]") {
	const auto cloud = createSphericalCloud(36, 18);
	savePointsPLY("sphere_36_18_cloud.ply", cloud);
	const auto mesh = measuredReconstruct(cloud, 0.3);
	CHECK(!mesh.empty());
	saveSTL("sphere_36_18_mesh.stl", cloud, mesh);
}

TEST_CASE("sphere_100_50", "[reconstruct]") {
	const auto cloud = createSphericalCloud(100, 50);
	savePointsPLY("sphere_100_50_cloud.ply", cloud);
	const auto mesh = measuredReconstruct(cloud, 0.1);
	CHECK(!mesh.empty());
	saveSTL("sphere_100_50_mesh.stl", cloud, mesh);
}

TEST_CASE("sphere_200_100", "[reconstruct]") {
	const auto cloud = createSphericalCloud(200, 100);
	savePointsPLY("sphere_200_100_cloud.ply", cloud);
	const auto mesh = measuredReconstruct(cloud, 0.04);
	CHECK(!mesh.empty());
	saveSTL("sphere_200_100_mesh.stl", cloud, mesh);
}

TEST_CASE("tetrahedron", "[reconstruct]") {
	const auto cloud = std::vector<Point>{{{0, 0, 0}, glm::normalize(glm::dvec3{-1, -1, -1})}, {{0, 1, 0}, glm::normalize(glm::dvec3{0, 1, 0})},
		{{1, 0, 0}, glm::normalize(glm::dvec3{1, 0, 0})}, {{0, 0, 1}, glm::normalize(glm::dvec3{0, 0, 1})}};
	savePointsPLY("tetrahedron_cloud.ply", cloud);
	const auto mesh = measuredReconstruct(cloud, 2);
	CHECK(!mesh.empty());
	saveSTL("tetrahedron_mesh.stl", cloud, mesh);
}

TEST_CASE("cube", "[reconstruct]") {
	const auto cloud = std::vector<Point>{
		{{-1, -1, -1}, glm::normalize(glm::dvec3{-1, -1, -1})},
		{{-1, +1, -1}, glm::normalize(glm::dvec3{-1, +1, -1})},
		{{+1, +1, -1}, glm::normalize(glm::dvec3{+1, +1, -1})},
		{{+1, -1, -1}, glm::normalize(glm::dvec3{+1, -1, -1})},
		{{-1, -1, +1}, glm::normalize(glm::dvec3{-1, -1, +1})},
		{{-1, +1, +1}, glm::normalize(glm::dvec3{-1, +1, +1})},
		{{+1, +1, +1}, glm::normalize(glm::dvec3{+1, +1, +1})},
		{{+1, -1, +1}, glm::normalize(glm::dvec3{+1, -1, +1})},
	};
	savePointsPLY("cube_cloud.ply", cloud);
	const auto mesh = measuredReconstruct(cloud, 2);
	CHECK(!mesh.empty());
	saveSTL("cube_mesh.stl", cloud, mesh);
}

TEST_CASE("bunny", "[reconstruct]") {
	const auto cloud = loadXYZ("../test/data/bunny.xyz");
	const auto mesh = measuredReconstruct(cloud, 0.002);
	CHECK(!mesh.empty());
	saveSTL("bunny_mesh.stl", cloud, mesh);
}