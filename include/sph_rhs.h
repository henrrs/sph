#pragma once

#include "neighbor_list.h"
#include "particle.h"
#include "vec3.h"

struct ParticleRHS {
	double densityDerivative = 0.0;
	Vec3 acceleration = { 0.0, 0.0, 0.0 };
	Vec3 xsphCorrection = { 0.0, 0.0, 0.0 };
};

ParticleRHS computeParticleRHS(int particleIndex,
	const ParticleSet& particles,
	const ParticleSet& wallParticles,
	const NeighborListCSR& neighbors);
