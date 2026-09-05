#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
#include <vector>

namespace bpa {
	struct Point {
		glm::dvec3 pos;
		glm::dvec3 normal;
	};

	// Three indices into the input point vector, counter-clockwise seen from outside.
	using Face = std::array<std::uint32_t, 3>;

	struct Options {
		// Only the nearest this many neighbours of a seed candidate are paired (0: all within
		// 2 x radius, the paper's unbounded search).
		std::size_t seedNeighbors = 100;
		// Drop connected components with fewer triangles than this (0: keep everything).
		std::size_t minComponent = 0;
	};

	auto reconstruct(const std::vector<Point>& points, double radius, const Options& options = {}) -> std::vector<Face>;

	// Several radii, in one pass each from the smallest (section 4.6 of the paper): the mesh of
	// one pass is kept, and its boundary edges whose triangles admit an empty ball of the next
	// radius resume pivoting with it, so that gaps the smaller ball could not cross are closed.
	auto reconstruct(const std::vector<Point>& points, std::vector<double> radii, const Options& options = {}) -> std::vector<Face>;
} // namespace bpa
