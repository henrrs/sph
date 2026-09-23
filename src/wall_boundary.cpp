#include "wall_boundary.h"

#include <cmath>

#include "equation_of_state.h"
#include "parameters.h"
#include "sph_kernels.h"
#include "vec3.h"

namespace {

constexpr double kSmallDistance = 1.0e-12;
const Vec3 kGravityVector = { 0.0, -sim::gravity, 0.0 };

[[nodiscard]] inline Vec3 particlePosition(const ParticleSet& particles, int index) noexcept {
	return { particles.posX[index], particles.posY[index], particles.posZ[index] };
}

[[nodiscard]] inline Vec3 particleVelocity(const ParticleSet& particles, int index) noexcept {
	return { particles.velX[index], particles.velY[index], particles.velZ[index] };
}

} // namespace

void resetWallParticleState(int wallParticleIndex, ParticleSet& wallParticles) {
	wallParticles.velX[wallParticleIndex] = 0.0;
	wallParticles.velY[wallParticleIndex] = 0.0;
	wallParticles.velZ[wallParticleIndex] = 0.0;
	wallParticles.previousVelX[wallParticleIndex] = 0.0;
	wallParticles.previousVelY[wallParticleIndex] = 0.0;
	wallParticles.previousVelZ[wallParticleIndex] = 0.0;
	wallParticles.pressure[wallParticleIndex] = 0.0;
	wallParticles.density[wallParticleIndex] = sim::densityOfFluid;
}

void updateWallParticleState(int wallParticleIndex,
	const ParticleSet& particles,
	ParticleSet& wallParticles,
	const NeighborListCSR& neighbors) {
	if (wallParticleIndex < 0 || wallParticleIndex + 1 >= static_cast<int>(neighbors.wallFluidOffsets.size())) {
		resetWallParticleState(wallParticleIndex, wallParticles);
		return;
	}

	const Vec3 xb = particlePosition(wallParticles, wallParticleIndex);
	Vec3 weightedVelocity = { 0.0, 0.0, 0.0 };
	double weightedPressure = 0.0;
	double hydrostaticCorrection = 0.0;
	double weightSum = 0.0;

	for (int k = neighbors.wallFluidOffsets[wallParticleIndex]; k < neighbors.wallFluidOffsets[wallParticleIndex + 1]; k++) {
		const int fluidIndex = neighbors.wallFluidIndices[k];
		const Vec3 xf = particlePosition(particles, fluidIndex);
		const Vec3 xbf = xb - xf;
		const double r2 = dot(xbf, xbf);
		const double wij = kernelFunction(std::sqrt(r2));
		const Vec3 vf = particleVelocity(particles, fluidIndex);

		weightSum += wij;
		weightedVelocity += wij * vf;
		weightedPressure += particles.pressure[fluidIndex] * wij;
		hydrostaticCorrection += particles.density[fluidIndex] * dot(kGravityVector, xbf) * wij;
	}

	if (weightSum <= kSmallDistance) {
		resetWallParticleState(wallParticleIndex, wallParticles);
		return;
	}

	const Vec3 wallVelocity = -1.0 * ((1.0 / weightSum) * weightedVelocity);

	wallParticles.velX[wallParticleIndex] = wallVelocity.x;
	wallParticles.velY[wallParticleIndex] = wallVelocity.y;
	wallParticles.velZ[wallParticleIndex] = wallVelocity.z;
	wallParticles.previousVelX[wallParticleIndex] = wallVelocity.x;
	wallParticles.previousVelY[wallParticleIndex] = wallVelocity.y;
	wallParticles.previousVelZ[wallParticleIndex] = wallVelocity.z;
	wallParticles.pressure[wallParticleIndex] = (weightedPressure + hydrostaticCorrection) / weightSum;
	wallParticles.density[wallParticleIndex] = calcDensityFromPressure(wallParticles.pressure[wallParticleIndex]);
}
