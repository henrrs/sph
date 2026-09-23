#include "simulation.h"

#include <algorithm>
#include <cstdio>
#include <iomanip>
#include <iostream>

#include "diagnostics.h"
#include "output_data.h"
#include "parameters.h"

SphSimulation::SphSimulation(SimulationConfig config)
	: config_(config),
	backend_(createSimulationBackend(config_)) {}

void SphSimulation::initialize() {
	backend_->initialize();
	if (config_.enableFrameDiagnostics) {
		writeHydrostaticInitializationDiagnostics(backend_->fluidParticles());
	}
}

int SphSimulation::stepsPerFrame() const noexcept {
	return std::max(1, config_.numberOfSteps / config_.numFrames);
}

double SphSimulation::physicalTimeAtStep(int step) const noexcept {
	return step * sim::dt;
}

double SphSimulation::elapsedSeconds(Clock::time_point start, Clock::time_point end) const {
	return std::chrono::duration<double>(end - start).count();
}

void SphSimulation::run() {
	const Clock::time_point programStart = Clock::now();
	const int frameStride = stepsPerFrame();
	int frameCounter = 0;

	for (int step = 0; step < config_.numberOfSteps; step++) {
		const Clock::time_point stepStart = Clock::now();
		double frameSeconds = 0.0;
		double frameDiagnosticsSeconds = 0.0;
		double framePreparationSeconds = 0.0;
		double frameWriteSeconds = 0.0;

		backend_->advanceStep();

		if (config_.printProgress) {
			std::printf("\n Iteration number: %d \n", step);
		}

		if ((step % frameStride) == 0) {
			const Clock::time_point frameStart = Clock::now();
			if (config_.printProgress) {
				std::printf("\n* Generating files. Number %d.\n", frameCounter);
			}

			const double physicalTime = physicalTimeAtStep(step);
			const BackendFrameDiagnostics frameDiagnostics = backend_->frameDiagnostics();
			if (config_.enableFrameDiagnostics) {
				const Clock::time_point diagnosticsStart = Clock::now();
				appendSimulationDiagnostics(backend_->fluidParticles(),
					step,
					physicalTime,
					frameDiagnostics.activeWallParticleCount,
					frameDiagnostics.fluidNeighborPairCount,
					frameDiagnostics.boundaryNeighborPairCount);
				frameDiagnosticsSeconds = elapsedSeconds(diagnosticsStart, Clock::now());
			}

			const Clock::time_point framePreparationStart = Clock::now();
			OutputFrameData frame = backend_->downloadFrameData();
			framePreparationSeconds = elapsedSeconds(framePreparationStart, Clock::now());

			const Clock::time_point frameWriteStart = Clock::now();
			outputWriter_.writeFrame(frame, "data", frameCounter, physicalTime);
			frameWriteSeconds = elapsedSeconds(frameWriteStart, Clock::now());

			frameCounter++;
			frameSeconds = elapsedSeconds(frameStart, Clock::now());
		}

		const double totalStepSeconds = elapsedSeconds(stepStart, Clock::now());
		if (config_.enableTimingDiagnostics) {
			const BackendStepTiming backendTiming = backend_->stepTiming();
			appendTimingDiagnostics(step,
				totalStepSeconds - frameSeconds,
				frameSeconds,
				totalStepSeconds);
			appendDetailedTimingDiagnostics(step,
				backendTiming,
				frameDiagnosticsSeconds,
				framePreparationSeconds,
				frameWriteSeconds,
				totalStepSeconds);
		}
	}

	const double timeTaken = elapsedSeconds(programStart, Clock::now());
	std::cout << "\nTime taken by program is : " << std::fixed << std::setprecision(5) << timeTaken;
	std::cout << " sec " << std::endl;
}
