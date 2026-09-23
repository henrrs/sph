#pragma once

#include <memory>

#include "output_data.h"
#include "particle.h"
#include "simulation_config.h"

struct BackendFrameDiagnostics {
	int activeWallParticleCount = 0;
	int fluidNeighborPairCount = 0;
	int boundaryNeighborPairCount = 0;
};

struct BackendStepTiming {
	double cellAssignmentSeconds = 0.0;
	double neighborBuildSeconds = 0.0;
	double wallActivationSeconds = 0.0;
	double pressureSeconds = 0.0;
	double wallUpdateSeconds = 0.0;
	double wallResetSeconds = 0.0;
	double rhsSeconds = 0.0;
	double stateCopySeconds = 0.0;
	double integrationSeconds = 0.0;
	double stateApplySeconds = 0.0;
	double cellClearSeconds = 0.0;
};

class SimulationBackend {
public:
	virtual ~SimulationBackend() = default;

	virtual void initialize() = 0;
	virtual void advanceStep() = 0;
	virtual OutputFrameData downloadFrameData() const = 0;
	virtual const ParticleSet& fluidParticles() const = 0;
	virtual const ParticleSet& wallParticles() const = 0;
	virtual BackendFrameDiagnostics frameDiagnostics() const = 0;
	virtual BackendStepTiming stepTiming() const = 0;
};

std::unique_ptr<SimulationBackend> createSimulationBackend(const SimulationConfig& config);
