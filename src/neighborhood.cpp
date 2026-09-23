#include "neighborhood.h"

#include <cmath>
#include <cstdio>
#include <vector>

#include "parameters.h"

namespace {

using namespace sim;
constexpr double khSquared = kh * kh;

inline double squaredDistance(double xi, double yi, double zi,
	double xj, double yj, double zj) noexcept {
	const double dx = xj - xi;
	const double dy = yj - yi;
	const double dz = zj - zi;
	return (dx * dx) + (dy * dy) + (dz * dz);
}

int coordinateToCellIndex(double coordinate, double domainMin, int cellCount) {
	const int index = static_cast<int>(std::floor((coordinate - domainMin) / kh));
	if (index < 0 || index >= cellCount) {
		const double domainMax = domainMin + (cellCount * kh);
		if (std::fabs(coordinate - domainMax) <= 1e-10) {
			return cellCount - 1;
		}
		return -1;
	}
	return index;
}

int computeCellTag(double x, double y, double z) {
	const int line = coordinateToCellIndex(x, domainMinX, numberOfCellsInX);
	const int column = coordinateToCellIndex(y, domainMinY, numberOfCellsInY);
	const int depth = coordinateToCellIndex(z, domainMinZ, numberOfCellsInZ);

	if (line < 0 || column < 0 || depth < 0) {
		return -1;
	}

	return depth + (line * numberOfCellsInZ) + (column * numberOfCellsInX * numberOfCellsInZ);
}

void countNeighborsForParticle(int particleIndex,
	const std::vector<Cell>& cells,
	const ParticleSet& particles,
	const ParticleSet& wallParticles,
	int& fluidCount,
	int& boundaryCount) {
	const double xi = particles.posX[particleIndex];
	const double yi = particles.posY[particleIndex];
	const double zi = particles.posZ[particleIndex];
	const int cellOfParticle = particles.currentCell[particleIndex];

	if (cellOfParticle < 0 || cellOfParticle >= static_cast<int>(cells.size())) {
		fluidCount = 0;
		boundaryCount = 0;
		return;
	}

	for (int neighborIdx = 0; neighborIdx < cells[cellOfParticle].neighborCount; neighborIdx++) {
		const int neighborCellTag = cells[cellOfParticle].neighborTags[neighborIdx];
		const Cell& neighborCell = cells[neighborCellTag];

		for (int fluidNeighborIndex : neighborCell.particles) {
			if (fluidNeighborIndex == particleIndex) {
				continue;
			}

			const double r2 = squaredDistance(xi, yi, zi,
				particles.posX[fluidNeighborIndex],
				particles.posY[fluidNeighborIndex],
				particles.posZ[fluidNeighborIndex]);
			if (r2 <= khSquared) {
				fluidCount++;
			}
		}

		for (int wallNeighborIndex : neighborCell.boundaryParticles) {
			const double r2 = squaredDistance(xi, yi, zi,
				wallParticles.posX[wallNeighborIndex],
				wallParticles.posY[wallNeighborIndex],
				wallParticles.posZ[wallNeighborIndex]);
			if (r2 <= khSquared) {
				boundaryCount++;
			}
		}
	}
}

void fillNeighborsForParticle(int particleIndex,
	const std::vector<Cell>& cells,
	const ParticleSet& particles,
	const ParticleSet& wallParticles,
	NeighborListCSR& neighbors) {
	const double xi = particles.posX[particleIndex];
	const double yi = particles.posY[particleIndex];
	const double zi = particles.posZ[particleIndex];
	const int cellOfParticle = particles.currentCell[particleIndex];

	if (cellOfParticle < 0 || cellOfParticle >= static_cast<int>(cells.size())) {
		return;
	}

	int fluidCursor = neighbors.fluidOffsets[particleIndex];
	int boundaryCursor = neighbors.boundaryOffsets[particleIndex];

	for (int neighborIdx = 0; neighborIdx < cells[cellOfParticle].neighborCount; neighborIdx++) {
		const int neighborCellTag = cells[cellOfParticle].neighborTags[neighborIdx];
		const Cell& neighborCell = cells[neighborCellTag];

		for (int fluidNeighborIndex : neighborCell.particles) {
			if (fluidNeighborIndex == particleIndex) {
				continue;
			}

			const double r2 = squaredDistance(xi, yi, zi,
				particles.posX[fluidNeighborIndex],
				particles.posY[fluidNeighborIndex],
				particles.posZ[fluidNeighborIndex]);
			if (r2 <= khSquared) {
				neighbors.fluidIndices[fluidCursor++] = fluidNeighborIndex;
			}
		}

		for (int wallNeighborIndex : neighborCell.boundaryParticles) {
			const double r2 = squaredDistance(xi, yi, zi,
				wallParticles.posX[wallNeighborIndex],
				wallParticles.posY[wallNeighborIndex],
				wallParticles.posZ[wallNeighborIndex]);
			if (r2 <= khSquared) {
				neighbors.boundaryIndices[boundaryCursor++] = wallNeighborIndex;
			}
		}
	}
}

} // namespace

void putWallParticlesInsideCells(std::vector<Cell>& cells, ParticleSet& particles) {
	for (std::size_t i = 0; i < particles.size(); i++) {
		const int cellTag = computeCellTag(particles.posX[i], particles.posY[i], particles.posZ[i]);
		particles.currentCell[i] = cellTag;
		if (cellTag >= 0) {
			cells[cellTag].boundaryParticles.push_back(static_cast<int>(i));
		}
	}
}

void putParticlesInsideCells(std::vector<Cell>& cells, ParticleSet& particles) {
	for (std::size_t i = 0; i < particles.size(); i++) {
		const int cellTag = computeCellTag(particles.posX[i], particles.posY[i], particles.posZ[i]);
		particles.currentCell[i] = cellTag;
		if (cellTag >= 0) {
			cells[cellTag].particles.push_back(static_cast<int>(i));
		}
	}
}

void buildNeighborsCSR(const std::vector<Cell>& cells,
	const ParticleSet& particles,
	const ParticleSet& wallParticles,
	NeighborListCSR& neighbors,
	NeighborBuildWorkspace& workspace) {
	const int particleCount = static_cast<int>(particles.size());
	const int wallParticleCount = static_cast<int>(wallParticles.size());
	workspace.fluidCounts.assign(particleCount, 0);
	workspace.boundaryCounts.assign(particleCount, 0);

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < particleCount; i++) {
		countNeighborsForParticle(i, cells, particles, wallParticles, workspace.fluidCounts[i], workspace.boundaryCounts[i]);
	}

	neighbors.resizeOffsets(particles.size());
	for (int i = 0; i < particleCount; i++) {
		neighbors.fluidOffsets[i + 1] = neighbors.fluidOffsets[i] + workspace.fluidCounts[i];
		neighbors.boundaryOffsets[i + 1] = neighbors.boundaryOffsets[i] + workspace.boundaryCounts[i];
	}

	neighbors.fluidIndices.resize(neighbors.fluidOffsets[particleCount]);
	neighbors.boundaryIndices.resize(neighbors.boundaryOffsets[particleCount]);

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < particleCount; i++) {
		fillNeighborsForParticle(i, cells, particles, wallParticles, neighbors);
	}

	workspace.wallFluidCounts.assign(wallParticleCount, 0);
	for (int wallIndex : neighbors.boundaryIndices) {
		workspace.wallFluidCounts[wallIndex]++;
	}

	neighbors.wallFluidOffsets.assign(wallParticleCount + 1, 0);
	neighbors.activeBoundaryIndices.clear();
	for (int i = 0; i < wallParticleCount; i++) {
		neighbors.wallFluidOffsets[i + 1] = neighbors.wallFluidOffsets[i] + workspace.wallFluidCounts[i];
		if (workspace.wallFluidCounts[i] > 0) {
			neighbors.activeBoundaryIndices.push_back(i);
		}
	}

	neighbors.wallFluidIndices.resize(neighbors.wallFluidOffsets[wallParticleCount]);
	workspace.wallFluidWriteOffsets = neighbors.wallFluidOffsets;
	for (int fluidIndex = 0; fluidIndex < particleCount; fluidIndex++) {
		for (int k = neighbors.boundaryOffsets[fluidIndex]; k < neighbors.boundaryOffsets[fluidIndex + 1]; k++) {
			const int wallIndex = neighbors.boundaryIndices[k];
			neighbors.wallFluidIndices[workspace.wallFluidWriteOffsets[wallIndex]++] = fluidIndex;
		}
	}
}

void takeParticlesOutOfCells(std::vector<Cell>& cells) {
	for (Cell& cell : cells) {
		cell.particles.clear();
	}
}

void clearParticleCellAssignments(ParticleSet& particles) {
	for (std::size_t i = 0; i < particles.size(); i++) {
		particles.currentCell[i] = -1;
	}
}
