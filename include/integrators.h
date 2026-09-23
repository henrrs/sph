#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include "neighbor_list.h"
#include "parameters.h"
#include "particle.h"
#include "vec3.h"

struct ParticleStateSoA {
	std::vector<double> posX;
	std::vector<double> posY;
	std::vector<double> posZ;
	std::vector<double> velX;
	std::vector<double> velY;
	std::vector<double> velZ;
	std::vector<double> previousVelX;
	std::vector<double> previousVelY;
	std::vector<double> previousVelZ;
	std::vector<double> accX;
	std::vector<double> accY;
	std::vector<double> accZ;
	std::vector<double> xsphX;
	std::vector<double> xsphY;
	std::vector<double> xsphZ;
	std::vector<double> density;
	std::vector<double> densityDerivation;
};

ParticleStateSoA createParticleStateBuffers(std::size_t count);
ParticleStateSoA snapshotParticleState(const ParticleSet& particles);
void copyParticleStateToBuffer(const ParticleSet& particles, ParticleStateSoA& state);
void integrateLeapfrogStep(int particleIndex, const ParticleStateSoA& state, ParticleStateSoA& nextState);
void integrateVelocityVerletMidpointPredictor(int particleIndex,
	const ParticleStateSoA& state,
	ParticleStateSoA& midpointState);
void integrateVelocityVerletFullPredictor(int particleIndex,
	const ParticleStateSoA& originalState,
	const ParticleStateSoA& midpointState,
	ParticleStateSoA& predictedState);
void integrateVelocityVerletCorrector(int particleIndex,
	const ParticleStateSoA& predictedState,
	const ParticleStateSoA& evaluatedState,
	ParticleStateSoA& correctedState);
void applyParticleState(const ParticleStateSoA& nextState, ParticleSet& particles);
