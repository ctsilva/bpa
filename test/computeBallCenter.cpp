#include <bpa.cpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace bpa;
using Catch::Matchers::WithinAbs;

namespace {
	void checkCenter(const dvec3& c, double x, double y, double z) {
		CHECK_THAT(c.x, WithinAbs(x, 1e-12));
		CHECK_THAT(c.y, WithinAbs(y, 1e-12));
		CHECK_THAT(c.z, WithinAbs(z, 1e-12));
	}
} // namespace

TEST_CASE("isosceles", "[computeBallCenter]") {
	MeshPoint a{dvec3{0, 0, 0}};
	MeshPoint b{dvec3{10, 0, 0}};
	MeshPoint c{dvec3{0, 10, 0}};
	MeshFace f{&a, &b, &c};
	const auto center = computeBallCenter(f, 10);
	REQUIRE(center);
	checkCenter(center.value(), 5, 5, std::sqrt(50.0));
}

TEST_CASE("isoscelesLargerRadius", "[computeBallCenter]") {
	MeshPoint a{dvec3{0, 0, 0}};
	MeshPoint b{dvec3{10, 0, 0}};
	MeshPoint c{dvec3{0, 10, 0}};
	MeshFace f{&a, &b, &c};
	const auto center = computeBallCenter(f, 100);
	REQUIRE(center);
	checkCenter(center.value(), 5, 5, std::sqrt(10000.0 - 50.0));
}

TEST_CASE("equilateral", "[computeBallCenter]") {
	MeshPoint a{dvec3{0, 0, 0}};
	MeshPoint b{dvec3{10, 0, 0}};
	MeshPoint c{dvec3{5, std::sqrt(3.0) * 5, 0}};
	MeshFace f{&a, &b, &c};
	const auto center = computeBallCenter(f, 10);
	REQUIRE(center);
	const auto circumradius = 10 / std::sqrt(3.0);
	checkCenter(center.value(), 5, 5 / std::sqrt(3.0), std::sqrt(100 - circumradius * circumradius));
}

TEST_CASE("radiusTooSmall", "[computeBallCenter]") {
	MeshPoint a{dvec3{0, 0, 0}};
	MeshPoint b{dvec3{10, 0, 0}};
	MeshPoint c{dvec3{0, 10, 0}};
	MeshFace f{&a, &b, &c};
	CHECK(!computeBallCenter(f, 1));
}

TEST_CASE("degenerate", "[computeBallCenter]") {
	MeshPoint a{dvec3{0, 0, 0}};
	MeshPoint b{dvec3{10, 0, 0}};
	MeshPoint c{dvec3{20, 0, 0}};
	MeshFace f{&a, &b, &c};
	CHECK(!computeBallCenter(f, 100));
}
