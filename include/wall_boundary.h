#pragma once

#include "neighbor_list.h"
#include "particle.h"

void resetWallParticleState(int wallParticleIndex, ParticleSet& wallParticles);
void updateWallParticleState(int wallParticleIndex,
	const ParticleSet& particles,
	ParticleSet& wallParticles,
	const NeighborListCSR& neighbors);
