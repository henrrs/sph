#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "equation_of_state.h"
#include "initialization.h"
#include "integrators.h"
#include "neighborhood.h"
#include "parameters.h"
#include "sph_kernels.h"
#include "sph_rhs.h"
#include "wall_boundary.h"

namespace {

struct TestContext {
	int failures = 0;

	void require(bool condition, const std::string& message) {
		if (condition) {
			std::cout << "[PASS] " << message << '\n';
			return;
		}
		failures++;
		std::cerr << "[FAIL] " << message << '\n';
	}

	void requireNear(double actual, double expected, double tolerance, const std::string& message) {
		const double error = std::abs(actual - expected);
		std::ostringstream details;
		details << message << " (actual=" << std::setprecision(15) << actual
			<< ", expected=" << expected << ", error=" << error
			<< ", tolerance=" << tolerance << ')';
		require(error <= tolerance, details.str());
	}
};

[[nodiscard]] bool nearlyEqual(double lhs, double rhs, double tolerance = 1.0e-12) {
	return std::abs(lhs - rhs) <= tolerance;
}

[[nodiscard]] int fluidIndex(int column, int line, int depth) {
	return column * sim::lineParticles * sim::depthParticles + line * sim::depthParticles + depth;
}

[[nodiscard]] Vec3 position(const ParticleSet& particles, int index) {
	return { particles.posX[index], particles.posY[index], particles.posZ[index] };
}

[[nodiscard]] Vec3 velocity(const ParticleSet& particles, int index) {
	return { particles.velX[index], particles.velY[index], particles.velZ[index] };
}

[[nodiscard]] double kernelPartition(int particleIndex,
	const ParticleSet& particles,
	const ParticleSet& wallParticles,
	const NeighborListCSR& neighbors) {
	double partition = (sim::particleMass / particles.density[particleIndex]) * kernelFunction(0.0);
	const Vec3 xi = position(particles, particleIndex);

	for (int k = neighbors.fluidOffsets[particleIndex]; k < neighbors.fluidOffsets[particleIndex + 1]; k++) {
		const int neighborIndex = neighbors.fluidIndices[k];
		const double distance = norm(xi - position(particles, neighborIndex));
		partition += (sim::particleMass / particles.density[neighborIndex]) * kernelFunction(distance);
	}

	for (int k = neighbors.boundaryOffsets[particleIndex]; k < neighbors.boundaryOffsets[particleIndex + 1]; k++) {
		const int wallIndex = neighbors.boundaryIndices[k];
		const double distance = norm(xi - position(wallParticles, wallIndex));
		partition += sim::referenceParticleVolume * kernelFunction(distance);
	}

	return partition;
}

void testWallGeometry(TestContext& test, const ParticleSet& wallParticles) {
	const int xSteps = static_cast<int>(std::lround(sim::xDomain / sim::space));
	const int ySteps = static_cast<int>(std::lround(sim::yDomain / sim::space));
	const int zSteps = static_cast<int>(std::lround(sim::zDomain / sim::space));
	const std::size_t expectedCount = static_cast<std::size_t>(sim::wallLayerCount) *
		((xSteps + (2 * sim::wallLayerCount)) * (zSteps + (2 * sim::wallLayerCount)) +
		 2 * (ySteps + sim::wallLayerCount) * (zSteps + (2 * sim::wallLayerCount)) +
		 2 * (ySteps + sim::wallLayerCount) * xSteps);

	test.require(wallParticles.size() == expectedCount, "wall particle count matches the disjoint face partition");

	bool hasBottom = false;
	bool hasLeft = false;
	bool hasRight = false;
	bool hasFront = false;
	bool hasBack = false;
	bool allCoordinatesAreHalfGrid = true;
	bool hasDuplicate = false;
	std::size_t duplicateIndex = 0;
	std::unordered_set<std::string> coordinates;
	coordinates.reserve(wallParticles.size());

	for (std::size_t i = 0; i < wallParticles.size(); i++) {
		const std::int64_t qx = std::llround(2.0 * wallParticles.posX[i] / sim::space);
		const std::int64_t qy = std::llround(2.0 * wallParticles.posY[i] / sim::space);
		const std::int64_t qz = std::llround(2.0 * wallParticles.posZ[i] / sim::space);
		allCoordinatesAreHalfGrid = allCoordinatesAreHalfGrid &&
			(std::abs(qx) % 2 == 1) && (std::abs(qy) % 2 == 1) && (std::abs(qz) % 2 == 1);

		const std::string key = std::to_string(qx) + ':' + std::to_string(qy) + ':' + std::to_string(qz);
		if (!coordinates.insert(key).second) {
			hasDuplicate = true;
			duplicateIndex = i;
			break;
		}

		hasBottom = hasBottom ||
			(nearlyEqual(wallParticles.posX[i], 0.5 * sim::space) &&
			 nearlyEqual(wallParticles.posY[i], -0.5 * sim::space) &&
			 nearlyEqual(wallParticles.posZ[i], 0.5 * sim::space));
		hasLeft = hasLeft ||
			(nearlyEqual(wallParticles.posX[i], -0.5 * sim::space) &&
			 nearlyEqual(wallParticles.posY[i], 0.5 * sim::space) &&
			 nearlyEqual(wallParticles.posZ[i], 0.5 * sim::space));
		hasRight = hasRight ||
			(nearlyEqual(wallParticles.posX[i], sim::xDomain + (0.5 * sim::space)) &&
			 nearlyEqual(wallParticles.posY[i], 0.5 * sim::space) &&
			 nearlyEqual(wallParticles.posZ[i], 0.5 * sim::space));
		hasFront = hasFront ||
			(nearlyEqual(wallParticles.posX[i], 0.5 * sim::space) &&
			 nearlyEqual(wallParticles.posY[i], 0.5 * sim::space) &&
			 nearlyEqual(wallParticles.posZ[i], -0.5 * sim::space));
		hasBack = hasBack ||
			(nearlyEqual(wallParticles.posX[i], 0.5 * sim::space) &&
			 nearlyEqual(wallParticles.posY[i], 0.5 * sim::space) &&
			 nearlyEqual(wallParticles.posZ[i], sim::zDomain + (0.5 * sim::space)));
	}

	test.require(!hasDuplicate,
		"wall face partitions do not overlap at edges or corners"
			+ (hasDuplicate ? " (first duplicate at index " + std::to_string(duplicateIndex) + ')' : ""));
	test.require(allCoordinatesAreHalfGrid, "all wall coordinates continue the fluid half-grid lattice");
	test.require(hasBottom && hasLeft && hasRight && hasFront && hasBack,
		"the five physical wall interfaces start one particle spacing from the adjacent fluid row");
}

void testKernelCoverageAndHydrostatics(TestContext& test,
	const ParticleSet& particles,
	const ParticleSet& wallParticles,
	const NeighborListCSR& neighbors) {
	const std::vector<std::pair<std::string, int>> samples = {
		{ "bottom", fluidIndex(0, sim::lineParticles / 2, sim::depthParticles / 2) },
		{ "bottom-front", fluidIndex(0, sim::lineParticles / 2, 0) },
		{ "bottom-left-front", fluidIndex(0, 0, 0) },
		{ "left", fluidIndex(sim::columnParticles / 2, 0, sim::depthParticles / 2) },
		{ "front", fluidIndex(sim::columnParticles / 2, sim::lineParticles / 2, 0) },
		{ "interior", fluidIndex(sim::columnParticles / 2, sim::lineParticles / 2, sim::depthParticles / 2) },
	};

	for (const auto& [name, index] : samples) {
		const double partition = kernelPartition(index, particles, wallParticles, neighbors);
		std::cout << "[METRIC] partition." << name << '=' << std::setprecision(12) << partition << '\n';
		test.require(std::abs(partition - 1.0) < 5.0e-3,
			"kernel partition is complete at sample " + name);

		const ParticleRHS rhs = computeParticleRHS(index, particles, wallParticles, neighbors);
		const double accelerationNorm = norm(rhs.acceleration);
		std::cout << "[METRIC] hydrostatic_acceleration." << name << '='
			<< std::setprecision(12) << accelerationNorm << " m/s^2\n";
		test.require(accelerationNorm < 0.5,
			"initial hydrostatic residual is small at sample " + name);
	}
}

void testFixedWallVolumeInContinuity(TestContext& test,
	ParticleSet particles,
	ParticleSet wallParticles,
	const NeighborListCSR& neighbors) {
	const int index = fluidIndex(0, sim::lineParticles / 2, sim::depthParticles / 2);
	particles.velX[index] = 0.25;
	particles.velY[index] = -1.0;
	particles.velZ[index] = 0.125;
	for (double& wallDensity : wallParticles.density) {
		wallDensity = 0.5 * sim::densityOfFluid;
	}

	const Vec3 xi = position(particles, index);
	const Vec3 vi = velocity(particles, index);
	const double rhoI = particles.density[index];
	double expectedFixedVolume = 0.0;
	double incorrectDynamicVolume = 0.0;

	for (int k = neighbors.fluidOffsets[index]; k < neighbors.fluidOffsets[index + 1]; k++) {
		const int neighborIndex = neighbors.fluidIndices[k];
		const Vec3 xij = xi - position(particles, neighborIndex);
		const Vec3 vij = vi - velocity(particles, neighborIndex);
		const Vec3 gradW = gradKernelFunction(norm(xij), xij);
		const double contribution =
			rhoI * (sim::particleMass / particles.density[neighborIndex]) * dot(vij, gradW);
		expectedFixedVolume += contribution;
		incorrectDynamicVolume += contribution;
	}

	for (int k = neighbors.boundaryOffsets[index]; k < neighbors.boundaryOffsets[index + 1]; k++) {
		const int wallIndex = neighbors.boundaryIndices[k];
		const Vec3 xij = xi - position(wallParticles, wallIndex);
		const Vec3 gradW = gradKernelFunction(norm(xij), xij);
		expectedFixedVolume += rhoI * sim::referenceParticleVolume * dot(vi, gradW);
		incorrectDynamicVolume +=
			rhoI * (sim::particleMass / wallParticles.density[wallIndex]) * dot(vi, gradW);
	}

	const ParticleRHS rhs = computeParticleRHS(index, particles, wallParticles, neighbors);
	const double tolerance = 1.0e-10 * std::max(1.0, std::abs(expectedFixedVolume));
	test.requireNear(rhs.densityDerivative, expectedFixedVolume, tolerance,
		"wall continuity uses the initial wall-particle volume");
	test.require(std::abs(rhs.densityDerivative - incorrectDynamicVolume) > 1.0,
		"wall continuity is independent of the perturbed hydrodynamic wall density");
}

void testAdamiPressurePair(TestContext& test) {
	ParticleSet particles;
	particles.resize(2);
	particles.posX = { 0.0, sim::space };
	particles.posY = { 0.0, 0.0 };
	particles.posZ = { 0.0, 0.0 };
	particles.velX = { 0.0, 0.0 };
	particles.velY = { 0.0, 0.0 };
	particles.velZ = { 0.0, 0.0 };
	particles.density = { 1010.0, 990.0 };
	particles.pressure = { 2000.0, -1000.0 };

	ParticleSet wallParticles;
	NeighborListCSR neighbors;
	neighbors.fluidOffsets = { 0, 1, 2 };
	neighbors.fluidIndices = { 1, 0 };
	neighbors.boundaryOffsets = { 0, 0, 0 };

	const Vec3 xij = position(particles, 0) - position(particles, 1);
	const Vec3 gradW = gradKernelFunction(norm(xij), xij);
	const double rhoI = particles.density[0];
	const double rhoJ = particles.density[1];
	const double pressureI = particles.pressure[0];
	const double pressureJ = particles.pressure[1];
	const double interparticlePressure = ((rhoJ * pressureI) + (rhoI * pressureJ)) / (rhoI + rhoJ);
	const double volumeI = sim::particleMass / rhoI;
	const double volumeJ = sim::particleMass / rhoJ;
	const double pairScale =
		((volumeI * volumeI) + (volumeJ * volumeJ)) * interparticlePressure / sim::particleMass;
	const Vec3 expectedAcceleration = (-pairScale * gradW) + Vec3 { 0.0, -sim::gravity, 0.0 };
	const ParticleRHS rhs = computeParticleRHS(0, particles, wallParticles, neighbors);

	const double tolerance = 1.0e-12 * std::max(1.0, norm(expectedAcceleration));
	test.requireNear(rhs.acceleration.x, expectedAcceleration.x, tolerance,
		"pressure force uses Adami's density-weighted interparticle pressure");
	test.requireNear(rhs.acceleration.y, expectedAcceleration.y, tolerance,
		"pressure-pair test preserves the body-force contribution");
	test.requireNear(rhs.acceleration.z, expectedAcceleration.z, tolerance,
		"pressure-pair test has no spurious transverse force");
}

void testVelocityVerletStages(TestContext& test) {
	ParticleStateSoA original = createParticleStateBuffers(1);
	ParticleStateSoA midpoint = createParticleStateBuffers(1);
	ParticleStateSoA predicted = createParticleStateBuffers(1);
	ParticleStateSoA evaluated = createParticleStateBuffers(1);
	ParticleStateSoA corrected = createParticleStateBuffers(1);

	original.posX[0] = 1.0;
	original.posY[0] = 2.0;
	original.posZ[0] = 3.0;
	original.velX[0] = 0.4;
	original.velY[0] = -0.2;
	original.velZ[0] = 0.1;
	original.accX[0] = 2.0;
	original.accY[0] = -4.0;
	original.accZ[0] = 6.0;
	original.xsphX[0] = 0.01;
	original.xsphY[0] = 0.02;
	original.xsphZ[0] = -0.03;
	original.density[0] = 1000.0;
	original.densityDerivation[0] = 8.0;

	integrateVelocityVerletMidpointPredictor(0, original, midpoint);
	const Vec3 expectedHalfVelocity = {
		original.velX[0] + (0.5 * sim::dt * original.accX[0]),
		original.velY[0] + (0.5 * sim::dt * original.accY[0]),
		original.velZ[0] + (0.5 * sim::dt * original.accZ[0])
	};
	test.requireNear(midpoint.velX[0], expectedHalfVelocity.x, 1.0e-15,
		"Velocity-Verlet predicts velocity at n+1/2");
	test.requireNear(midpoint.posX[0],
		original.posX[0] + (0.5 * sim::dt * (expectedHalfVelocity.x + original.xsphX[0])),
		1.0e-15,
		"Velocity-Verlet transports the first half-step with the initial XSPH field");
	test.requireNear(midpoint.density[0],
		original.density[0] + (0.5 * sim::dt * original.densityDerivation[0]),
		1.0e-12,
		"Velocity-Verlet predicts density at n+1/2");

	midpoint.xsphX[0] = -0.02;
	midpoint.xsphY[0] = 0.04;
	midpoint.xsphZ[0] = 0.01;
	midpoint.densityDerivation[0] = -6.0;
	integrateVelocityVerletFullPredictor(0, original, midpoint, predicted);
	test.requireNear(predicted.posX[0],
		midpoint.posX[0] + (0.5 * sim::dt * (midpoint.velX[0] + midpoint.xsphX[0])),
		1.0e-15,
		"Velocity-Verlet transports the second half-step with midpoint XSPH");
	test.requireNear(predicted.density[0],
		original.density[0] + (sim::dt * midpoint.densityDerivation[0]),
		1.0e-12,
		"Velocity-Verlet integrates density with the midpoint derivative");

	evaluated.accX[0] = -3.0;
	evaluated.accY[0] = 5.0;
	evaluated.accZ[0] = -7.0;
	evaluated.densityDerivation[0] = 11.0;
	evaluated.xsphX[0] = 0.5;
	evaluated.xsphY[0] = 0.25;
	evaluated.xsphZ[0] = -0.125;
	integrateVelocityVerletCorrector(0, predicted, evaluated, corrected);
	test.requireNear(corrected.velX[0],
		predicted.velX[0] + (0.5 * sim::dt * evaluated.accX[0]),
		1.0e-15,
		"Velocity-Verlet corrects velocity with acceleration at n+1");
	test.requireNear(corrected.density[0], predicted.density[0], 1.0e-15,
		"Velocity-Verlet corrector preserves midpoint-integrated density");
	test.requireNear(corrected.densityDerivation[0], evaluated.densityDerivation[0], 1.0e-15,
		"Velocity-Verlet caches the evaluated RHS for the next step");
}

} // namespace

int main() {
	TestContext test;

	std::vector<Cell> cells;
	ParticleSet particles;
	ParticleSet wallParticles;
	NeighborListCSR neighbors;
	NeighborBuildWorkspace workspace;

	startCells(cells);
	startParticles(particles);
	startWalls(wallParticles);
	testWallGeometry(test, wallParticles);

	putWallParticlesInsideCells(cells, wallParticles);
	putParticlesInsideCells(cells, particles);
	buildNeighborsCSR(cells, particles, wallParticles, neighbors, workspace);
	updatePressureField(particles.density, particles.pressure);
	for (int wallIndex : neighbors.activeBoundaryIndices) {
		updateWallParticleState(wallIndex, particles, wallParticles, neighbors);
	}

	testKernelCoverageAndHydrostatics(test, particles, wallParticles, neighbors);
	testFixedWallVolumeInContinuity(test, particles, wallParticles, neighbors);
	testAdamiPressurePair(test);
	testVelocityVerletStages(test);

	if (test.failures != 0) {
		std::cerr << "\n" << test.failures << " physics consistency test(s) failed.\n";
		return 1;
	}

	std::cout << "\nAll physics consistency tests passed.\n";
	return 0;
}
