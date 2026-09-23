#include "sph_rhs.h"

#include <cmath>

#include "parameters.h"
#include "sph_kernels.h"

namespace {

[[nodiscard]] inline double square(double value) noexcept {
	return value * value;
}

struct PairInteraction {
	Vec3 xij = { 0.0, 0.0, 0.0 };
	Vec3 vij = { 0.0, 0.0, 0.0 };
	double r2 = 0.0;
	double r = 0.0;
	Vec3 gradW = { 0.0, 0.0, 0.0 };
};

[[nodiscard]] inline Vec3 particlePosition(const ParticleSet& particles, int index) noexcept {
	return { particles.posX[index], particles.posY[index], particles.posZ[index] };
}

[[nodiscard]] inline Vec3 particleVelocity(const ParticleSet& particles, int index) noexcept {
	return { particles.velX[index], particles.velY[index], particles.velZ[index] };
}

[[nodiscard]] inline double pressureGradientPair(double pressureI, double rhoI, double pressureJ, double rhoJ) noexcept {
	const double interparticlePressure = ((rhoJ * pressureI) + (rhoI * pressureJ)) / (rhoI + rhoJ);
	const double volumeI = sim::particleMass / rhoI;
	const double volumeJ = sim::particleMass / rhoJ;
	return ((square(volumeI) + square(volumeJ)) / sim::particleMass) * interparticlePressure;
}

[[nodiscard]] PairInteraction makePairInteraction(const Vec3& xi,
	const Vec3& xj,
	const Vec3& vi,
	const Vec3& vj) {
	PairInteraction pair;
	pair.xij = xi - xj;
	pair.vij = vi - vj;
	pair.r2 = dot(pair.xij, pair.xij);
	pair.r = std::sqrt(pair.r2);
	pair.gradW = gradKernelFunction(pair.r, pair.xij);
	return pair;
}

[[nodiscard]] double calcArtificialViscosityPair(const PairInteraction& pair,
	double rhoI,
	double rhoJ) noexcept {
	const double scalarProduct = dot(pair.vij, pair.xij);
	if (scalarProduct < 0.0) {
		const double phi_ij = (scalarProduct * sim::h) / (pair.r2 + (0.01 * square(sim::h)));
		const double rhoij = (rhoI + rhoJ) / 2.0;
		return ((sim::alpha * phi_ij * sim::speedOfSound) + (sim::beta * square(phi_ij))) / rhoij;
	}
	return 0.0;
}

void addDensityContribution(double& densityDerivative,
	const PairInteraction& pair,
	double rhoI,
	double neighborVolume) {
	densityDerivative += rhoI * neighborVolume * dot(pair.vij, pair.gradW);
}

void addAccelerationContribution(Vec3& acceleration,
	const PairInteraction& pair,
	double pressureI,
	double rhoI,
	double pressureJ,
	double rhoJ,
	bool includeArtificialViscosity) {
	const double artificialViscosity = includeArtificialViscosity ? calcArtificialViscosityPair(pair, rhoI, rhoJ) : 0.0;
	const double pressureGradient = pressureGradientPair(pressureI, rhoI, pressureJ, rhoJ);
	const double pairScale = pressureGradient + (artificialViscosity * sim::particleMass);
	acceleration -= pairScale * pair.gradW;
}

} // namespace

ParticleRHS computeParticleRHS(int particleIndex,
	const ParticleSet& particles,
	const ParticleSet& wallParticles,
	const NeighborListCSR& neighbors) {
	ParticleRHS rhs;
	const Vec3 xi = particlePosition(particles, particleIndex);
	const Vec3 vi = particleVelocity(particles, particleIndex);
	const Vec3 staticWallVelocity = { 0.0, 0.0, 0.0 };
	const double pressureI = particles.pressure[particleIndex];
	const double rhoI = particles.density[particleIndex];

	for (int k = neighbors.fluidOffsets[particleIndex]; k < neighbors.fluidOffsets[particleIndex + 1]; k++) {
		const int neighborIndex = neighbors.fluidIndices[k];
		const Vec3 xj = particlePosition(particles, neighborIndex);
		const Vec3 vj = particleVelocity(particles, neighborIndex);
		const PairInteraction pair = makePairInteraction(xi, xj, vi, vj);
		const double wij = kernelFunction(pair.r);
		const double rhoJ = particles.density[neighborIndex];

		addDensityContribution(rhs.densityDerivative, pair, rhoI, sim::particleMass / rhoJ);
		addAccelerationContribution(rhs.acceleration, pair, pressureI, rhoI, particles.pressure[neighborIndex], rhoJ, true);

		const double rhoij = 0.5 * (rhoI + rhoJ);
		rhs.xsphCorrection += (sim::particleMass * wij / rhoij) * pair.vij;
	}

	for (int k = neighbors.boundaryOffsets[particleIndex]; k < neighbors.boundaryOffsets[particleIndex + 1]; k++) {
		const int wallIndex = neighbors.boundaryIndices[k];
		const Vec3 xj = particlePosition(wallParticles, wallIndex);
		const Vec3 vj = particleVelocity(wallParticles, wallIndex);
		const PairInteraction pair = makePairInteraction(xi, xj, vi, vj);

		const PairInteraction continuityPair = { pair.xij, vi - staticWallVelocity, pair.r2, pair.r, pair.gradW };
		addDensityContribution(rhs.densityDerivative, continuityPair, rhoI, sim::referenceParticleVolume);
		// Adami et al. use artificial viscosity only for fluid-fluid pairs.
		addAccelerationContribution(rhs.acceleration,
			pair,
			pressureI,
			rhoI,
			wallParticles.pressure[wallIndex],
			wallParticles.density[wallIndex],
			false);
	}

	rhs.xsphCorrection = sim::xsphConstant * rhs.xsphCorrection;
	rhs.acceleration.y -= sim::gravity;
	return rhs;
}
