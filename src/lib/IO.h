#pragma once

#include "bpa.h"

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

namespace bpa {
	// Point cloud readers. The format is chosen by the extension: ".xyz" holds "x y z nx ny nz"
	// per line; ".noff" (or ".off" with normals) is the Object File Format with a "NOFF" header,
	// "nv nf ne", then "x y z nx ny nz" per vertex; faces, if any, are ignored.
	inline auto loadXYZ(const std::filesystem::path& path) -> std::vector<Point> {
		std::ifstream f{path};
		if (!f)
			throw std::runtime_error("Failed to read file " + path.string());
		std::vector<Point> result;
		Point p;
		while (f >> p.pos.x >> p.pos.y >> p.pos.z >> p.normal.x >> p.normal.y >> p.normal.z)
			result.push_back(p);
		return result;
	}

	inline auto loadNOFF(const std::filesystem::path& path) -> std::vector<Point> {
		std::ifstream f{path};
		if (!f)
			throw std::runtime_error("Failed to read file " + path.string());
		std::string magic;
		std::size_t nv = 0, nf = 0, ne = 0;
		f >> magic >> nv >> nf >> ne;
		if (magic != "NOFF")
			throw std::runtime_error(path.string() + ": expected a NOFF header (vertices with normals)");
		std::vector<Point> result(nv);
		for (auto& p : result)
			f >> p.pos.x >> p.pos.y >> p.pos.z >> p.normal.x >> p.normal.y >> p.normal.z;
		if (!f)
			throw std::runtime_error(path.string() + ": truncated vertex list");
		return result;
	}

	inline auto loadPoints(const std::filesystem::path& path) -> std::vector<Point> {
		const auto ext = path.extension().string();
		if (ext == ".noff" || ext == ".off")
			return loadNOFF(path);
		return loadXYZ(path);
	}

	// Mesh writers. OFF keeps the input vertices in their order, so the faces index them
	// directly; STL is a binary triangle soup with float32 positions.
	inline void saveOFF(const std::filesystem::path& path, const std::vector<Point>& points, const std::vector<Face>& faces) {
		if (path.has_parent_path())
			create_directories(path.parent_path());
		std::ofstream f{path};
		f.precision(9);
		f << "OFF\n" << points.size() << ' ' << faces.size() << " 0\n";
		for (const auto& p : points)
			f << p.pos.x << ' ' << p.pos.y << ' ' << p.pos.z << '\n';
		for (const auto& t : faces)
			f << "3 " << t[0] << ' ' << t[1] << ' ' << t[2] << '\n';
	}

	inline void saveSTL(const std::filesystem::path& path, const std::vector<Point>& points, const std::vector<Face>& faces) {
		if (path.has_parent_path())
			create_directories(path.parent_path());
		std::ofstream f{path, std::ios::binary};
		const char header[80] = "binary STL written by bpa";
		f.write(header, sizeof(header));
		const auto count = static_cast<std::uint32_t>(faces.size());
		f.write(reinterpret_cast<const char*>(&count), sizeof(count));
		const std::uint16_t attributeCount = 0;
		for (const auto& t : faces) {
			const glm::vec3 v[3] = {points[t[0]].pos, points[t[1]].pos, points[t[2]].pos};
			const auto normal = glm::normalize(glm::cross(v[1] - v[0], v[2] - v[0]));
			f.write(reinterpret_cast<const char*>(&normal), sizeof(normal));
			f.write(reinterpret_cast<const char*>(v), sizeof(v));
			f.write(reinterpret_cast<const char*>(&attributeCount), sizeof(attributeCount));
		}
	}

	inline void saveMesh(const std::filesystem::path& path, const std::vector<Point>& points, const std::vector<Face>& faces) {
		if (path.extension() == ".stl")
			saveSTL(path, points, faces);
		else
			saveOFF(path, points, faces);
	}

	inline void savePointsPLY(const std::filesystem::path& path, const std::vector<Point>& points) {
		if (path.has_parent_path())
			create_directories(path.parent_path());
		std::ofstream f{path, std::ios::binary};
		f << "ply\nformat binary_little_endian 1.0\nelement vertex " << points.size()
		  << "\nproperty float x\nproperty float y\nproperty float z\nproperty float nx\nproperty float ny\nproperty float nz\nend_header\n";
		f.write(reinterpret_cast<const char*>(points.data()), points.size() * sizeof(points[0]));
	}
} // namespace bpa
