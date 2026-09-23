#include "initialization.h"

#include <cmath>

namespace {

using namespace sim;

int cellIndexFromCoordinates(int line, int column, int depth) {
	return depth + (line * numberOfCellsInZ) + (column * numberOfCellsInX * numberOfCellsInZ);
}

void appendNeighbor(Cell& cell, int neighborTag) {
	cell.neighborTags[cell.neighborCount] = neighborTag;
	cell.neighborCount++;
}

void findNeighborsOfACell(std::vector<Cell>& cells, int cellTag) {
	Cell& cell = cells[cellTag];
	cell.neighborCount = 0;

	for (int deltaColumn = -1; deltaColumn <= 1; deltaColumn++) {
		const int neighborColumn = cell.column + deltaColumn;
		if (neighborColumn < 0 || neighborColumn >= numberOfCellsInY) {
			continue;
		}

		for (int deltaLine = -1; deltaLine <= 1; deltaLine++) {
			const int neighborLine = cell.line + deltaLine;
			if (neighborLine < 0 || neighborLine >= numberOfCellsInX) {
				continue;
			}

			for (int deltaDepth = -1; deltaDepth <= 1; deltaDepth++) {
				const int neighborDepth = cell.depth + deltaDepth;
				if (neighborDepth < 0 || neighborDepth >= numberOfCellsInZ) {
					continue;
				}

				appendNeighbor(cell, cellIndexFromCoordinates(neighborLine, neighborColumn, neighborDepth));
			}
		}
	}
}

void appendParticle(ParticleSet& particles,
	double x,
	double y,
	double z,
	int particleType,
	int particleTag,
	double density,
	double pressure) {
	particles.posX.push_back(x);
	particles.posY.push_back(y);
	particles.posZ.push_back(z);
	particles.velX.push_back(0.0);
	particles.velY.push_back(0.0);
	particles.velZ.push_back(0.0);
	particles.previousVelX.push_back(0.0);
	particles.previousVelY.push_back(0.0);
	particles.previousVelZ.push_back(0.0);
	particles.accX.push_back(0.0);
	particles.accY.push_back(0.0);
	particles.accZ.push_back(0.0);
	particles.xsphX.push_back(0.0);
	particles.xsphY.push_back(0.0);
	particles.xsphZ.push_back(0.0);
	particles.density.push_back(density);
	particles.densityDerivation.push_back(0.0);
	particles.pressure.push_back(pressure);
	particles.type.push_back(particleType);
	particles.tag.push_back(particleTag);
	particles.currentCell.push_back(-1);
}

} // namespace

void startCells(std::vector<Cell>& cells) {
	using namespace sim;

	cells.clear();
	cells.reserve(totalNumCells);

	int cellTag = 0;
	for (int column = 0; column < numberOfCellsInY; column++) {
		const double startCellInY = domainMinY + (column * kh);
		const double endCellInY = startCellInY + kh;

		for (int line = 0; line < numberOfCellsInX; line++) {
			const double startCellInX = domainMinX + (line * kh);
			const double endCellInX = startCellInX + kh;

			for (int depth = 0; depth < numberOfCellsInZ; depth++) {
				const double startCellInZ = domainMinZ + (depth * kh);
				const double endCellInZ = startCellInZ + kh;

				Cell cell;
				cell.x0 = startCellInX;
				cell.x1 = endCellInX;
				cell.y0 = startCellInY;
				cell.y1 = endCellInY;
				cell.z0 = startCellInZ;
				cell.z1 = endCellInZ;
				cell.tag = cellTag;
				cell.line = line;
				cell.column = column;
				cell.depth = depth;
				cell.particles.reserve(128);
				cell.boundaryParticles.reserve(128);
				cells.push_back(cell);
				findNeighborsOfACell(cells, cellTag);
				cellTag++;
			}
		}
	}
}

void startParticles(ParticleSet& particles) {
	using namespace sim;

	particles.clear();
	particles.resize(numFluidParticles);

	int tag = 0;
	for (int j = 0; j < columnParticles; j++) {
		const double y = (0.5 * space) + (j * space);
		for (int i = 0; i < lineParticles; i++) {
			const double x = (0.5 * space) + (i * space);
			for (int k = 0; k < depthParticles; k++) {
				const double z = (0.5 * space) + (k * space);
				const double density = calcInitialDensity(y);

				particles.posX[tag] = x;
				particles.posY[tag] = y;
				particles.posZ[tag] = z;
				particles.velX[tag] = 0.0;
				particles.velY[tag] = 0.0;
				particles.velZ[tag] = 0.0;
				particles.previousVelX[tag] = 0.0;
				particles.previousVelY[tag] = 0.0;
				particles.previousVelZ[tag] = 0.0;
				particles.accX[tag] = 0.0;
				particles.accY[tag] = 0.0;
				particles.accZ[tag] = 0.0;
				particles.xsphX[tag] = 0.0;
				particles.xsphY[tag] = 0.0;
				particles.xsphZ[tag] = 0.0;
				particles.density[tag] = density;
				particles.densityDerivation[tag] = 0.0;
				particles.pressure[tag] = 0.0;
				particles.type[tag] = fluidParticleType;
				particles.tag[tag] = tag;
				particles.currentCell[tag] = -1;
				tag++;
			}
		}
	}
}

void startWalls(ParticleSet& wallParticles) {
	using namespace sim;

	wallParticles.clear();

	const int xSteps = static_cast<int>(std::lround(xDomain / space));
	const int ySteps = static_cast<int>(std::lround(yDomain / space));
	const int zSteps = static_cast<int>(std::lround(zDomain / space));
	const std::size_t estimatedWallParticles = static_cast<std::size_t>(wallLayerCount) *
		((xSteps + (2 * wallLayerCount)) * (zSteps + (2 * wallLayerCount)) +
		 2 * (ySteps + wallLayerCount) * (zSteps + (2 * wallLayerCount)) +
		 2 * (ySteps + wallLayerCount) * xSteps);
	wallParticles.reserve(estimatedWallParticles);

	int tag = numFluidParticles;
	for (int layer = 1; layer <= wallLayerCount; layer++) {
		const double y = -(static_cast<double>(layer) - 0.5) * space;
		for (int ix = -wallLayerCount; ix < xSteps + wallLayerCount; ix++) {
			const double x = (static_cast<double>(ix) + 0.5) * space;
			for (int iz = -wallLayerCount; iz < zSteps + wallLayerCount; iz++) {
				const double z = (static_cast<double>(iz) + 0.5) * space;
				appendParticle(wallParticles, x, y, z, wallParticleType, tag, densityOfFluid, 0.0);
				tag++;
			}
		}
	}

	for (int layer = 1; layer <= wallLayerCount; layer++) {
		const double normalOffset = (static_cast<double>(layer) - 0.5) * space;
		const double leftX = -normalOffset;
		const double rightX = xDomain + normalOffset;
		for (int iy = 0; iy < ySteps + wallLayerCount; iy++) {
			const double y = (static_cast<double>(iy) + 0.5) * space;
			for (int iz = -wallLayerCount; iz < zSteps + wallLayerCount; iz++) {
				const double z = (static_cast<double>(iz) + 0.5) * space;
				appendParticle(wallParticles, leftX, y, z, wallParticleType, tag, densityOfFluid, 0.0);
				tag++;
				appendParticle(wallParticles, rightX, y, z, wallParticleType, tag, densityOfFluid, 0.0);
				tag++;
			}
		}
	}

	for (int layer = 1; layer <= wallLayerCount; layer++) {
		const double normalOffset = (static_cast<double>(layer) - 0.5) * space;
		const double frontZ = -normalOffset;
		const double backZ = zDomain + normalOffset;
		for (int iy = 0; iy < ySteps + wallLayerCount; iy++) {
			const double y = (static_cast<double>(iy) + 0.5) * space;
			for (int ix = 0; ix < xSteps; ix++) {
				const double x = (static_cast<double>(ix) + 0.5) * space;
				appendParticle(wallParticles, x, y, frontZ, wallParticleType, tag, densityOfFluid, 0.0);
				tag++;
				appendParticle(wallParticles, x, y, backZ, wallParticleType, tag, densityOfFluid, 0.0);
				tag++;
			}
		}
	}
}

double calcInitialDensity(double y) {
	const double B = ((sim::speedOfSound * sim::speedOfSound) * sim::densityOfFluid) / sim::gamma;
	return sim::densityOfFluid * std::pow(1.0 + (sim::densityOfFluid * sim::gravity * (sim::fluidHeight - y)) / B, 1.0 / sim::gamma);
}
