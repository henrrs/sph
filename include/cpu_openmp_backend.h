#pragma once

#include <vector>

#include "cell.h"
#include "integrators.h"
#include "neighbor_list.h"
#include "simulation_backend.h"

class CpuOpenMpBackend final : public SimulationBackend {
public:
	explicit CpuOpenMpBackend(const SimulationConfig& config);

	void initialize() override;
	void advanceStep() override;
	OutputFrameData downloadFrameData() const override;
	const ParticleSet& fluidParticles() const override;
	const ParticleSet& wallParticles() const override;
	BackendFrameDiagnostics frameDiagnostics() const override;
	BackendStepTiming stepTiming() const override;

private:
	void refreshActiveWallParticles();
	void rebuildFluidNeighborhood();
	void updateParticleDynamics();
	void clearFluidCellState();
	void resetStepTiming();

	SimulationConfig config_;
	ParticleSet fluidParticles_;
	ParticleSet wallParticles_;
	std::vector<Cell> cells_;
	NeighborListCSR neighbors_;
	NeighborBuildWorkspace neighborWorkspace_;
	std::vector<unsigned char> activeWallParticleMarkers_;
	std::vector<int> activeWallParticleIndices_;
	std::vector<int> previousActiveWallParticleIndices_;
	ParticleStateSoA state_;
	ParticleStateSoA midpointState_;
	ParticleStateSoA nextState_;
	bool dynamicsInitialized_ = false;
	BackendStepTiming stepTiming_;
};
