#pragma once

#include <vector>

#include "cell.h"
#include "neighbor_list.h"
#include "particle.h"

void putParticlesInsideCells(std::vector<Cell>& cells, ParticleSet& particles);
void putWallParticlesInsideCells(std::vector<Cell>& cells, ParticleSet& particles);
void buildNeighborsCSR(const std::vector<Cell>& cells,
	const ParticleSet& particles,
	const ParticleSet& wallParticles,
	NeighborListCSR& neighbors,
	NeighborBuildWorkspace& workspace);
void takeParticlesOutOfCells(std::vector<Cell>& cells);
void clearParticleCellAssignments(ParticleSet& particles);
