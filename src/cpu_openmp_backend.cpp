#include "cpu_openmp_backend.h"

#include <algorithm>
#include <chrono>

#include <omp.h>

#include "equation_of_state.h"
#include "initialization.h"
#include "neighborhood.h"
#include "output_data.h"
#include "parameters.h"
#include "sph_rhs.h"
#include "wall_boundary.h"

namespace {

using Clock = std::chrono::steady_clock;

double elapsedSeconds(Clock::time_point start, Clock::time_point end) {
	return std::chrono::duration<double>(end - start).count();
}

} // namespace

CpuOpenMpBackend::CpuOpenMpBackend(const SimulationConfig& config)
	: config_(config) {}

void CpuOpenMpBackend::initialize() {
	startCells(cells_);
	startParticles(fluidParticles_);
	startWalls(wallParticles_);
	putWallParticlesInsideCells(cells_, wallParticles_);
	activeWallParticleMarkers_.resize(wallParticles_.size(), 0);
	state_ = createParticleStateBuffers(fluidParticles_.size());
	midpointState_ = createParticleStateBuffers(fluidParticles_.size());
	nextState_ = createParticleStateBuffers(fluidParticles_.size());
	dynamicsInitialized_ = false;
}

void CpuOpenMpBackend::refreshActiveWallParticles() {
	for (int wallIndex : previousActiveWallParticleIndices_) {
		activeWallParticleMarkers_[wallIndex] = 0;
	}
	activeWallParticleIndices_ = neighbors_.activeBoundaryIndices;
	for (int wallIndex : activeWallParticleIndices_) {
		activeWallParticleMarkers_[wallIndex] = 1;
	}
}

void CpuOpenMpBackend::rebuildFluidNeighborhood() {
	const Clock::time_point assignmentStart = Clock::now();
	putParticlesInsideCells(cells_, fluidParticles_);
	stepTiming_.cellAssignmentSeconds += elapsedSeconds(assignmentStart, Clock::now());

	const Clock::time_point neighborStart = Clock::now();
	buildNeighborsCSR(cells_, fluidParticles_, wallParticles_, neighbors_, neighborWorkspace_);
	stepTiming_.neighborBuildSeconds += elapsedSeconds(neighborStart, Clock::now());
}

void CpuOpenMpBackend::updateParticleDynamics() {
	const Clock::time_point activationStart = Clock::now();
	refreshActiveWallParticles();
	stepTiming_.wallActivationSeconds += elapsedSeconds(activationStart, Clock::now());

	Clock::time_point phaseStart = Clock::now();

	#pragma omp parallel
	{
		#pragma omp for schedule(static)
		for (int i = 0; i < static_cast<int>(fluidParticles_.size()); i++) {
			fluidParticles_.pressure[i] = calcPressure(fluidParticles_.density[i]);
		}

		#pragma omp single
		{
			stepTiming_.pressureSeconds += elapsedSeconds(phaseStart, Clock::now());
			phaseStart = Clock::now();
		}

		#pragma omp for schedule(guided, 64)
		for (int i = 0; i < static_cast<int>(activeWallParticleIndices_.size()); i++) {
			updateWallParticleState(activeWallParticleIndices_[i], fluidParticles_, wallParticles_, neighbors_);
		}

		#pragma omp single
		{
			stepTiming_.wallUpdateSeconds += elapsedSeconds(phaseStart, Clock::now());
			phaseStart = Clock::now();
		}

		#pragma omp for schedule(static)
		for (int i = 0; i < static_cast<int>(previousActiveWallParticleIndices_.size()); i++) {
			const int wallIndex = previousActiveWallParticleIndices_[i];
			if (activeWallParticleMarkers_[wallIndex] == 0) {
				resetWallParticleState(wallIndex, wallParticles_);
			}
		}

		#pragma omp single
		{
			stepTiming_.wallResetSeconds += elapsedSeconds(phaseStart, Clock::now());
			phaseStart = Clock::now();
		}

		#pragma omp for schedule(guided, 64)
		for (int i = 0; i < static_cast<int>(fluidParticles_.size()); i++) {
			const ParticleRHS rhs = computeParticleRHS(i, fluidParticles_, wallParticles_, neighbors_);
			fluidParticles_.densityDerivation[i] = rhs.densityDerivative;
			fluidParticles_.accX[i] = rhs.acceleration.x;
			fluidParticles_.accY[i] = rhs.acceleration.y;
			fluidParticles_.accZ[i] = rhs.acceleration.z;
			fluidParticles_.xsphX[i] = rhs.xsphCorrection.x;
			fluidParticles_.xsphY[i] = rhs.xsphCorrection.y;
			fluidParticles_.xsphZ[i] = rhs.xsphCorrection.z;
		}

		#pragma omp single
		{
			stepTiming_.rhsSeconds += elapsedSeconds(phaseStart, Clock::now());
		}
	}

	previousActiveWallParticleIndices_ = activeWallParticleIndices_;
}

void CpuOpenMpBackend::clearFluidCellState() {
	const Clock::time_point clearStart = Clock::now();
	takeParticlesOutOfCells(cells_);
	clearParticleCellAssignments(fluidParticles_);
	stepTiming_.cellClearSeconds += elapsedSeconds(clearStart, Clock::now());
}

void CpuOpenMpBackend::resetStepTiming() {
	stepTiming_ = {};
}

void CpuOpenMpBackend::advanceStep() {
	resetStepTiming();

	if constexpr (sim::useVelocityVerletIntegrator) {
		if (!dynamicsInitialized_) {
			rebuildFluidNeighborhood();
			updateParticleDynamics();
			dynamicsInitialized_ = true;
		}

		Clock::time_point phaseStart = Clock::now();
		copyParticleStateToBuffer(fluidParticles_, state_);
		stepTiming_.stateCopySeconds += elapsedSeconds(phaseStart, Clock::now());

		phaseStart = Clock::now();
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < static_cast<int>(fluidParticles_.size()); i++) {
			integrateVelocityVerletMidpointPredictor(i, state_, nextState_);
		}
		stepTiming_.integrationSeconds += elapsedSeconds(phaseStart, Clock::now());

		phaseStart = Clock::now();
		applyParticleState(nextState_, fluidParticles_);
		stepTiming_.stateApplySeconds += elapsedSeconds(phaseStart, Clock::now());
		clearFluidCellState();

		rebuildFluidNeighborhood();
		updateParticleDynamics();

		phaseStart = Clock::now();
		copyParticleStateToBuffer(fluidParticles_, midpointState_);
		stepTiming_.stateCopySeconds += elapsedSeconds(phaseStart, Clock::now());

		phaseStart = Clock::now();
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < static_cast<int>(fluidParticles_.size()); i++) {
			integrateVelocityVerletFullPredictor(i, state_, midpointState_, nextState_);
		}
		stepTiming_.integrationSeconds += elapsedSeconds(phaseStart, Clock::now());

		phaseStart = Clock::now();
		applyParticleState(nextState_, fluidParticles_);
		stepTiming_.stateApplySeconds += elapsedSeconds(phaseStart, Clock::now());
		clearFluidCellState();

		rebuildFluidNeighborhood();
		updateParticleDynamics();

		phaseStart = Clock::now();
		copyParticleStateToBuffer(fluidParticles_, midpointState_);
		stepTiming_.stateCopySeconds += elapsedSeconds(phaseStart, Clock::now());

		phaseStart = Clock::now();
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < static_cast<int>(fluidParticles_.size()); i++) {
			integrateVelocityVerletCorrector(i, nextState_, midpointState_, state_);
		}
		stepTiming_.integrationSeconds += elapsedSeconds(phaseStart, Clock::now());

		phaseStart = Clock::now();
		applyParticleState(state_, fluidParticles_);
		stepTiming_.stateApplySeconds += elapsedSeconds(phaseStart, Clock::now());
		clearFluidCellState();
	}
	else {
		rebuildFluidNeighborhood();
		updateParticleDynamics();

		Clock::time_point phaseStart = Clock::now();
		copyParticleStateToBuffer(fluidParticles_, state_);
		stepTiming_.stateCopySeconds += elapsedSeconds(phaseStart, Clock::now());

		phaseStart = Clock::now();
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < static_cast<int>(fluidParticles_.size()); i++) {
			integrateLeapfrogStep(i, state_, nextState_);
		}
		stepTiming_.integrationSeconds += elapsedSeconds(phaseStart, Clock::now());

		phaseStart = Clock::now();
		applyParticleState(nextState_, fluidParticles_);
		stepTiming_.stateApplySeconds += elapsedSeconds(phaseStart, Clock::now());
		clearFluidCellState();
	}
}

OutputFrameData CpuOpenMpBackend::downloadFrameData() const {
	return collectOutputFrameData(fluidParticles_, wallParticles_);
}

const ParticleSet& CpuOpenMpBackend::fluidParticles() const {
	return fluidParticles_;
}

const ParticleSet& CpuOpenMpBackend::wallParticles() const {
	return wallParticles_;
}

BackendFrameDiagnostics CpuOpenMpBackend::frameDiagnostics() const {
	return {
		static_cast<int>(activeWallParticleIndices_.size()),
		static_cast<int>(neighbors_.fluidIndices.size()),
		static_cast<int>(neighbors_.boundaryIndices.size()),
	};
}

BackendStepTiming CpuOpenMpBackend::stepTiming() const {
	return stepTiming_;
}
