#include <IO.h>
#include <bpa.h>

#include <chrono>
#include <iostream>

int main(int argc, char* argv[]) try {
	if (argc != 3 && argc != 4) {
		std::cerr << "Usage: " << argv[0] << " <points.xyz|points.noff> <radius> [<mesh.off|mesh.stl>]\n";
		return 1;
	}

	const std::filesystem::path inputFile = argv[1];
	const double radius = std::stod(argv[2]);
	const std::filesystem::path outputFile = argc == 4 ? argv[3] : inputFile.string() + ".off";

	const auto points = bpa::loadPoints(inputFile);
	const auto t0 = std::chrono::steady_clock::now();
	const auto faces = bpa::reconstruct(points, radius);
	const auto seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
	bpa::saveMesh(outputFile, points, faces);
	std::cout << points.size() << " points, " << faces.size() << " triangles, time: " << seconds << " s\n";

	return 0;
} catch (const std::exception& e) {
	std::cerr << "Error: " << e.what() << '\n';
	return 2;
}
