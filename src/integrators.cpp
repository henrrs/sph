#include "integrators.h"

#include <algorithm>
#include <span>

namespace {

using namespace sim;

void copyField(std::span<const double> source, std::span<double> destination) {
	std::copy(source.begin(), source.end(), destination.begin());
}

} // namespace

ParticleStateSoA createParticleStateBuffers(std::size_t count) {
	ParticleStateSoA state;
	state.posX.resize(count);
	state.posY.resize(count);
	state.posZ.resize(count);
	state.velX.resize(count);
	state.velY.resize(count);
	state.velZ.resize(count);
	state.previousVelX.resize(count);
	state.previousVelY.resize(count);
	state.previousVelZ.resize(count);
	state.accX.resize(count);
	state.accY.resize(count);
	state.accZ.resize(count);
	state.xsphX.resize(count);
	state.xsphY.resize(count);
	state.xsphZ.resize(count);
	state.density.resize(count);
	state.densityDerivation.resize(count);
	return state;
}

ParticleStateSoA snapshotParticleState(const ParticleSet& particles) {
	ParticleStateSoA state = createParticleStateBuffers(particles.size());
	copyParticleStateToBuffer(particles, state);
	return state;
}

void copyParticleStateToBuffer(const ParticleSet& particles, ParticleStateSoA& state) {
	if (state.posX.size() != particles.size()) {
		state = createParticleStateBuffers(particles.size());
	}
	copyField(particles.posX, state.posX);
	copyField(particles.posY, state.posY);
	copyField(particles.posZ, state.posZ);
	copyField(particles.velX, state.velX);
	copyField(particles.velY, state.velY);
	copyField(particles.velZ, state.velZ);
	copyField(particles.previousVelX, state.previousVelX);
	copyField(particles.previousVelY, state.previousVelY);
	copyField(particles.previousVelZ, state.previousVelZ);
	copyField(particles.accX, state.accX);
	copyField(particles.accY, state.accY);
	copyField(particles.accZ, state.accZ);
	copyField(particles.xsphX, state.xsphX);
	copyField(particles.xsphY, state.xsphY);
	copyField(particles.xsphZ, state.xsphZ);
	copyField(particles.density, state.density);
	copyField(particles.densityDerivation, state.densityDerivation);
}

void integrateLeapfrogStep(int particleIndex, const ParticleStateSoA& state, ParticleStateSoA& nextState) {
	const Vec3 correction = { state.xsphX[particleIndex], state.xsphY[particleIndex], state.xsphZ[particleIndex] };

	const double newPreviousVelocityX = state.previousVelX[particleIndex] + (dt * state.accX[particleIndex]);
	const double newPreviousVelocityY = state.previousVelY[particleIndex] + (dt * state.accY[particleIndex]);
	const double newPreviousVelocityZ = state.previousVelZ[particleIndex] + (dt * state.accZ[particleIndex]);

	const double newVelocityX = 0.5 * (newPreviousVelocityX + state.previousVelX[particleIndex]);
	const double newVelocityY = 0.5 * (newPreviousVelocityY + state.previousVelY[particleIndex]);
	const double newVelocityZ = 0.5 * (newPreviousVelocityZ + state.previousVelZ[particleIndex]);

	nextState.previousVelX[particleIndex] = newPreviousVelocityX;
	nextState.previousVelY[particleIndex] = newPreviousVelocityY;
	nextState.previousVelZ[particleIndex] = newPreviousVelocityZ;
	nextState.velX[particleIndex] = newVelocityX;
	nextState.velY[particleIndex] = newVelocityY;
	nextState.velZ[particleIndex] = newVelocityZ;
	nextState.posX[particleIndex] = state.posX[particleIndex] + (dt * (newPreviousVelocityX + correction.x));
	nextState.posY[particleIndex] = state.posY[particleIndex] + (dt * (newPreviousVelocityY + correction.y));
	nextState.posZ[particleIndex] = state.posZ[particleIndex] + (dt * (newPreviousVelocityZ + correction.z));
	nextState.density[particleIndex] = state.density[particleIndex] + (dt * state.densityDerivation[particleIndex]);
	nextState.densityDerivation[particleIndex] = state.densityDerivation[particleIndex];
	nextState.accX[particleIndex] = state.accX[particleIndex];
	nextState.accY[particleIndex] = state.accY[particleIndex];
	nextState.accZ[particleIndex] = state.accZ[particleIndex];
}

void integrateVelocityVerletMidpointPredictor(int particleIndex,
	const ParticleStateSoA& state,
	ParticleStateSoA& midpointState) {
	const Vec3 correction = { state.xsphX[particleIndex], state.xsphY[particleIndex], state.xsphZ[particleIndex] };

	const Vec3 halfVelocity = {
		state.velX[particleIndex] + (0.5 * dt * state.accX[particleIndex]),
		state.velY[particleIndex] + (0.5 * dt * state.accY[particleIndex]),
		state.velZ[particleIndex] + (0.5 * dt * state.accZ[particleIndex])
	};

	midpointState.previousVelX[particleIndex] = halfVelocity.x;
	midpointState.previousVelY[particleIndex] = halfVelocity.y;
	midpointState.previousVelZ[particleIndex] = halfVelocity.z;
	midpointState.velX[particleIndex] = halfVelocity.x;
	midpointState.velY[particleIndex] = halfVelocity.y;
	midpointState.velZ[particleIndex] = halfVelocity.z;
	midpointState.posX[particleIndex] = state.posX[particleIndex] + (0.5 * dt * (halfVelocity.x + correction.x));
	midpointState.posY[particleIndex] = state.posY[particleIndex] + (0.5 * dt * (halfVelocity.y + correction.y));
	midpointState.posZ[particleIndex] = state.posZ[particleIndex] + (0.5 * dt * (halfVelocity.z + correction.z));
	midpointState.density[particleIndex] = state.density[particleIndex] + (0.5 * dt * state.densityDerivation[particleIndex]);
	midpointState.densityDerivation[particleIndex] = state.densityDerivation[particleIndex];
	midpointState.accX[particleIndex] = state.accX[particleIndex];
	midpointState.accY[particleIndex] = state.accY[particleIndex];
	midpointState.accZ[particleIndex] = state.accZ[particleIndex];
	midpointState.xsphX[particleIndex] = state.xsphX[particleIndex];
	midpointState.xsphY[particleIndex] = state.xsphY[particleIndex];
	midpointState.xsphZ[particleIndex] = state.xsphZ[particleIndex];
}

void integrateVelocityVerletFullPredictor(int particleIndex,
	const ParticleStateSoA& originalState,
	const ParticleStateSoA& midpointState,
	ParticleStateSoA& predictedState) {
	const Vec3 midpointCorrection = {
		midpointState.xsphX[particleIndex],
		midpointState.xsphY[particleIndex],
		midpointState.xsphZ[particleIndex]
	};

	predictedState.previousVelX[particleIndex] = midpointState.velX[particleIndex];
	predictedState.previousVelY[particleIndex] = midpointState.velY[particleIndex];
	predictedState.previousVelZ[particleIndex] = midpointState.velZ[particleIndex];
	predictedState.velX[particleIndex] = midpointState.velX[particleIndex];
	predictedState.velY[particleIndex] = midpointState.velY[particleIndex];
	predictedState.velZ[particleIndex] = midpointState.velZ[particleIndex];
	predictedState.posX[particleIndex] =
		midpointState.posX[particleIndex] + (0.5 * dt * (midpointState.velX[particleIndex] + midpointCorrection.x));
	predictedState.posY[particleIndex] =
		midpointState.posY[particleIndex] + (0.5 * dt * (midpointState.velY[particleIndex] + midpointCorrection.y));
	predictedState.posZ[particleIndex] =
		midpointState.posZ[particleIndex] + (0.5 * dt * (midpointState.velZ[particleIndex] + midpointCorrection.z));
	predictedState.density[particleIndex] =
		originalState.density[particleIndex] + (dt * midpointState.densityDerivation[particleIndex]);
	predictedState.densityDerivation[particleIndex] = midpointState.densityDerivation[particleIndex];
	predictedState.accX[particleIndex] = midpointState.accX[particleIndex];
	predictedState.accY[particleIndex] = midpointState.accY[particleIndex];
	predictedState.accZ[particleIndex] = midpointState.accZ[particleIndex];
	predictedState.xsphX[particleIndex] = midpointState.xsphX[particleIndex];
	predictedState.xsphY[particleIndex] = midpointState.xsphY[particleIndex];
	predictedState.xsphZ[particleIndex] = midpointState.xsphZ[particleIndex];
}

void integrateVelocityVerletCorrector(int particleIndex,
	const ParticleStateSoA& predictedState,
	const ParticleStateSoA& evaluatedState,
	ParticleStateSoA& correctedState) {
	const Vec3 correctedVelocity = {
		predictedState.velX[particleIndex] + (0.5 * dt * evaluatedState.accX[particleIndex]),
		predictedState.velY[particleIndex] + (0.5 * dt * evaluatedState.accY[particleIndex]),
		predictedState.velZ[particleIndex] + (0.5 * dt * evaluatedState.accZ[particleIndex])
	};

	correctedState.velX[particleIndex] = correctedVelocity.x;
	correctedState.velY[particleIndex] = correctedVelocity.y;
	correctedState.velZ[particleIndex] = correctedVelocity.z;
	correctedState.previousVelX[particleIndex] = correctedVelocity.x;
	correctedState.previousVelY[particleIndex] = correctedVelocity.y;
	correctedState.previousVelZ[particleIndex] = correctedVelocity.z;
	correctedState.posX[particleIndex] = predictedState.posX[particleIndex];
	correctedState.posY[particleIndex] = predictedState.posY[particleIndex];
	correctedState.posZ[particleIndex] = predictedState.posZ[particleIndex];
	correctedState.density[particleIndex] = predictedState.density[particleIndex];
	correctedState.densityDerivation[particleIndex] = evaluatedState.densityDerivation[particleIndex];
	correctedState.accX[particleIndex] = evaluatedState.accX[particleIndex];
	correctedState.accY[particleIndex] = evaluatedState.accY[particleIndex];
	correctedState.accZ[particleIndex] = evaluatedState.accZ[particleIndex];
	correctedState.xsphX[particleIndex] = evaluatedState.xsphX[particleIndex];
	correctedState.xsphY[particleIndex] = evaluatedState.xsphY[particleIndex];
	correctedState.xsphZ[particleIndex] = evaluatedState.xsphZ[particleIndex];
}

void applyParticleState(const ParticleStateSoA& nextState, ParticleSet& particles) {
	copyField(nextState.posX, particles.posX);
	copyField(nextState.posY, particles.posY);
	copyField(nextState.posZ, particles.posZ);
	copyField(nextState.velX, particles.velX);
	copyField(nextState.velY, particles.velY);
	copyField(nextState.velZ, particles.velZ);
	copyField(nextState.previousVelX, particles.previousVelX);
	copyField(nextState.previousVelY, particles.previousVelY);
	copyField(nextState.previousVelZ, particles.previousVelZ);
	copyField(nextState.density, particles.density);
}
