#include <IO.h>
#include <bpa.h>
#include <chrono>
#include <cstring>
#include <iostream>
#include <string>

namespace {
	const char* usage = R"(usage: bpa <points> <radius> [<mesh>] [options]

  <points>   "x y z nx ny nz" per line (.xyz), or a NOFF file (.noff, .off)
  <radius>   the ball radius, in the units of the points, or a comma-separated list of
             radii for one pass each, smallest first
  <mesh>     .off (vertices in input order; the default, next to the input) or .stl

  --min-component N   drop connected components of fewer than N triangles
  --seed-neighbors N  pair a seed candidate only with its N nearest neighbours
                      (default 100; 0 pairs every point within 2 x radius)
)";
}

int main(int argc, char* argv[]) try {
	std::vector<std::string> positional;
	bpa::Options options;
	for (auto i = 1; i < argc; i++) {
		const std::string arg = argv[i];
		if (arg == "--min-component" && i + 1 < argc)
			options.minComponent = std::stoul(argv[++i]);
		else if (arg == "--seed-neighbors" && i + 1 < argc)
			options.seedNeighbors = std::stoul(argv[++i]);
		else if (arg == "-h" || arg == "--help" || arg.starts_with("-")) {
			std::cerr << usage;
			return arg == "-h" || arg == "--help" ? 0 : 1;
		} else
			positional.push_back(arg);
	}
	if (positional.size() != 2 && positional.size() != 3) {
		std::cerr << usage;
		return 1;
	}

	const std::filesystem::path inputFile = positional[0];
	std::vector<double> radii;
	for (std::size_t from = 0; from <= positional[1].size();) {
		const auto to = std::min(positional[1].find(',', from), positional[1].size());
		radii.push_back(std::stod(positional[1].substr(from, to - from)));
		from = to + 1;
	}
	const std::filesystem::path outputFile = positional.size() == 3 ? positional[2] : inputFile.string() + ".off";

	const auto points = bpa::loadPoints(inputFile);
	const auto t0 = std::chrono::steady_clock::now();
	const auto faces = bpa::reconstruct(points, radii, options);
	const auto seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
	bpa::saveMesh(outputFile, points, faces);
	std::cout << points.size() << " points, " << faces.size() << " triangles, time: " << seconds << " s\n";

	return 0;
} catch (const std::exception& e) {
	std::cerr << "Error: " << e.what() << '\n';
	return 2;
}
