#pragma once

#include <chrono>
#include <memory>

#include "output_writer.h"
#include "simulation_backend.h"
#include "simulation_config.h"

class SphSimulation {
public:
	explicit SphSimulation(SimulationConfig config);

	void initialize();
	void run();

private:
	using Clock = std::chrono::steady_clock;

	[[nodiscard]] int stepsPerFrame() const noexcept;
	[[nodiscard]] double physicalTimeAtStep(int step) const noexcept;
	[[nodiscard]] double elapsedSeconds(Clock::time_point start, Clock::time_point end) const;

	SimulationConfig config_;
	std::unique_ptr<SimulationBackend> backend_;
	SimulationOutputWriter outputWriter_;
};
