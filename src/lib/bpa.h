#pragma once

#include <array>
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

	auto reconstruct(const std::vector<Point>& points, double radius) -> std::vector<Face>;
} // namespace bpa
