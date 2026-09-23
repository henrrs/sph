#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <span>
#include <sstream>

#include "legacy_vtk_writer.h"
#include "output_data.h"
#include "parameters.h"

namespace {

namespace fs = std::filesystem;
using namespace sim;

fs::path buildOutputPath(const fs::path& subdirectory, const std::string& name, int width, int frameNumber, const char* extension) {
	fs::create_directories(subdirectory);
	std::stringstream ss;
	ss << name << "_" << std::setfill('0') << std::setw(width) << frameNumber << extension;
	return subdirectory / ss.str();
}

std::ofstream openOutputFile(const fs::path& path) {
	std::ofstream outFile(path, std::ios::out);
	if (!outFile) {
		std::cerr << "File couldn't be opened: " << path << '\n';
		std::exit(1);
	}
	return outFile;
}

void writePointBlock(std::ofstream& outFile,
	std::span<const double> posX,
	std::span<const double> posY,
	std::span<const double> posZ,
	int precision) {
	outFile << std::fixed << std::setprecision(precision);
	for (std::size_t i = 0; i < posX.size(); i++) {
		outFile << posX[i] << "  " << posY[i] << "  " << posZ[i] << '\n';
	}
}

void writePointBlock(std::ofstream& outFile,
	const std::vector<double>& points,
	int precision) {
	outFile << std::fixed << std::setprecision(precision);
	for (std::size_t i = 0; i < points.size(); i += 3) {
		outFile << points[i] << "  " << points[i + 1] << "  " << points[i + 2] << '\n';
	}
}

void writeScalarBlock(std::ofstream& outFile, const std::vector<double>& values) {
	for (double value : values) {
		outFile << value << '\n';
	}
}

void writeVectorBlock(std::ofstream& outFile, const std::vector<double>& velocity) {
	for (std::size_t i = 0; i < velocity.size(); i += 3) {
		outFile << velocity[i] << "  " << velocity[i + 1] << "  " << velocity[i + 2] << '\n';
	}
}

} // namespace

void createFile(const ParticleSet& particles, std::string name, int frameNumber) {
	const int precision = 5;
	const int width = 5;
	std::ofstream outFile = openOutputFile(buildOutputPath("output/mat", name, width, frameNumber, ".mat"));

	outFile << particles.size() << '\n';
	outFile << std::fixed << std::setprecision(precision);
	for (std::size_t i = 0; i < particles.size(); i++) {
		outFile << 0.0 << " " << 0.0 << " " << 0.0 << '\n';
	}
}

void createFilePositionMomentum(const ParticleSet& particles, std::string name, int frameNumber) {
	const int precision = 5;
	const int width = 5;
	const ParticleScalarsView view = makeParticleScalarsView(particles);
	std::ofstream outFile = openOutputFile(buildOutputPath("output/mat", name, width, frameNumber, ".mat"));

	outFile << 2 * particles.size() << '\n';
	writePointBlock(outFile, view.posX, view.posY, view.posZ, precision);

	outFile << std::fixed << std::setprecision(precision);
	for (std::size_t i = 0; i < particles.size(); i++) {
		outFile << particleMass * view.posX[i] << " " << particleMass * view.posY[i] << " " << particleMass * view.posZ[i] << '\n';
	}
}

void createPositionFile(const ParticleSet& particles, std::string name, int frameNumber) {
	const int precision = 5;
	const int width = 5;
	const ParticleScalarsView view = makeParticleScalarsView(particles);
	std::ofstream outFile = openOutputFile(buildOutputPath("output/mat", name, width, frameNumber, ".mat"));
	writePointBlock(outFile, view.posX, view.posY, view.posZ, precision);
}

void createMomentumFile(const ParticleSet& particles, std::string name, int frameNumber) {
	const int precision = 5;
	const int width = 5;
	const ParticleScalarsView view = makeParticleScalarsView(particles);
	std::ofstream outFile = openOutputFile(buildOutputPath("output/mat", name, width, frameNumber, ".mat"));

	outFile << particles.size() << '\n';
	outFile << std::fixed << std::setprecision(precision);
	for (std::size_t i = 0; i < particles.size(); i++) {
		outFile << particleMass * view.posX[i] << " " << particleMass * view.posY[i] << " " << particleMass * view.posZ[i] << '\n';
	}
}

void exportMeshLegacyVTK(const OutputFrameData& frame, std::string name, int frameNumber) {
	const int precision = 10;
	const int width = 2;
	const std::size_t totalPoints = frame.pointCount();
	std::ofstream outFile = openOutputFile(buildOutputPath("output/vtk", name, width, frameNumber, ".vtk"));

	outFile << "# vtk DataFile Version 2.0 \nUnstructured Grid \nASCII \nDATASET UNSTRUCTURED_GRID\n";
	outFile << "POINTS " << totalPoints << " double\n";
	writePointBlock(outFile, frame.points, precision);

	outFile << "CELLS " << totalPoints << "  " << 2 * totalPoints << '\n';
	for (std::size_t i = 0; i < totalPoints; i++) {
		outFile << 1 << "  " << i << '\n';
	}

	outFile << "CELL_TYPES " << totalPoints << '\n';
	for (std::size_t i = 0; i < totalPoints; i++) {
		outFile << 1 << '\n';
	}

	outFile << "POINT_DATA " << totalPoints << "\nSCALARS pressure double \nLOOKUP_TABLE default\n";
	writeScalarBlock(outFile, frame.pressure);

	outFile << "VECTORS velocity double\n";
	writeVectorBlock(outFile, frame.velocity);
}
