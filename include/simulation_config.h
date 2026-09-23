#pragma once

#include "parameters.h"

enum class BackendKind {
	CpuOpenMp,
};

struct SimulationConfig {
	int numberOfSteps = sim::numberOfSteps;
	int numFrames = sim::numFrames;
	BackendKind backend = BackendKind::CpuOpenMp;
	bool enableFrameDiagnostics = true;
	bool enableTimingDiagnostics = true;
	bool printProgress = true;
};

SimulationConfig makeDefaultSimulationConfig();
SimulationConfig parseSimulationConfig(int argc, char** argv);
