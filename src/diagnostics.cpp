#include "diagnostics.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>

#include "parameters.h"
#include "equation_of_state.h"

namespace {

namespace fs = std::filesystem;
using namespace sim;

fs::path diagnosticsDirectory() {
	fs::path directory = "output/diagnostics";
	fs::create_directories(directory);
	return directory;
}

std::ofstream openDiagnosticsFile(const fs::path& path, std::ios_base::openmode mode) {
	std::ofstream outFile(path, mode);
	if (!outFile) {
		throw std::runtime_error("Could not open diagnostics file: " + path.string());
	}
	return outFile;
}

struct FieldStats {
	double minValue = std::numeric_limits<double>::max();
	double maxValue = std::numeric_limits<double>::lowest();

	void add(double value) noexcept {
		minValue = std::min(minValue, value);
		maxValue = std::max(maxValue, value);
	}
};

} // namespace

void writeHydrostaticInitializationDiagnostics(const ParticleSet& particles) {
	const fs::path path = diagnosticsDirectory() / "hydrostatic_initial.csv";
	std::ofstream outFile = openDiagnosticsFile(path, std::ios::out | std::ios::trunc);

	outFile << "tag,y,density,eos_pressure,expected_pressure,pressure_error\n";
	outFile << std::fixed << std::setprecision(12);
	for (std::size_t i = 0; i < particles.size(); i++) {
		const double expectedPressure = densityOfFluid * gravity * (fluidHeight - particles.posY[i]);
		const double eosPressure = calcPressure(particles.density[i]);
		outFile << particles.tag[i] << ','
			<< particles.posY[i] << ','
			<< particles.density[i] << ','
			<< eosPressure << ','
			<< expectedPressure << ','
			<< (eosPressure - expectedPressure) << '\n';
	}
}

void appendSimulationDiagnostics(const ParticleSet& particles,
	int step,
	double physicalTime,
	int activeWallParticleCount,
	int fluidNeighborPairCount,
	int boundaryNeighborPairCount) {
	const fs::path path = diagnosticsDirectory() / "simulation.csv";
	const bool shouldWriteHeader = !fs::exists(path);
	std::ofstream outFile = openDiagnosticsFile(path, std::ios::out | std::ios::app);

	if (shouldWriteHeader) {
		outFile << "step,time,min_x,max_x,min_y,max_y,min_z,max_z,"
			<< "min_density,max_density,min_pressure,max_pressure,max_speed,"
			<< "below_bottom,above_top,outside_x,outside_z,"
			<< "active_wall_particles,fluid_neighbor_pairs,boundary_neighbor_pairs\n";
	}

	FieldStats xStats;
	FieldStats yStats;
	FieldStats zStats;
	FieldStats densityStats;
	FieldStats pressureStats;
	double maxSpeed = 0.0;
	int belowBottom = 0;
	int aboveTop = 0;
	int outsideX = 0;
	int outsideZ = 0;

	for (std::size_t i = 0; i < particles.size(); i++) {
		const double x = particles.posX[i];
		const double y = particles.posY[i];
		const double z = particles.posZ[i];
		const double speed = std::sqrt((particles.velX[i] * particles.velX[i])
			+ (particles.velY[i] * particles.velY[i])
			+ (particles.velZ[i] * particles.velZ[i]));

		xStats.add(x);
		yStats.add(y);
		zStats.add(z);
		densityStats.add(particles.density[i]);
		pressureStats.add(particles.pressure[i]);
		maxSpeed = std::max(maxSpeed, speed);

		belowBottom += y < 0.0 ? 1 : 0;
		aboveTop += y > yDomain ? 1 : 0;
		outsideX += (x < 0.0 || x > xDomain) ? 1 : 0;
		outsideZ += (z < 0.0 || z > zDomain) ? 1 : 0;
	}

	outFile << std::fixed << std::setprecision(12)
		<< step << ','
		<< physicalTime << ','
		<< xStats.minValue << ','
		<< xStats.maxValue << ','
		<< yStats.minValue << ','
		<< yStats.maxValue << ','
		<< zStats.minValue << ','
		<< zStats.maxValue << ','
		<< densityStats.minValue << ','
		<< densityStats.maxValue << ','
		<< pressureStats.minValue << ','
		<< pressureStats.maxValue << ','
		<< maxSpeed << ','
		<< belowBottom << ','
		<< aboveTop << ','
		<< outsideX << ','
		<< outsideZ << ','
		<< activeWallParticleCount << ','
		<< fluidNeighborPairCount << ','
		<< boundaryNeighborPairCount << '\n';
}

void appendTimingDiagnostics(int step,
	double computeSeconds,
	double frameSeconds,
	double totalSeconds) {
	const fs::path path = diagnosticsDirectory() / "timing.csv";
	const bool shouldWriteHeader = !fs::exists(path);
	std::ofstream outFile = openDiagnosticsFile(path, std::ios::out | std::ios::app);

	if (shouldWriteHeader) {
		outFile << "step,compute_seconds,frame_seconds,total_seconds\n";
	}

	outFile << std::fixed << std::setprecision(12)
		<< step << ','
		<< computeSeconds << ','
		<< frameSeconds << ','
		<< totalSeconds << '\n';
}

void appendDetailedTimingDiagnostics(int step,
	const BackendStepTiming& backendTiming,
	double frameDiagnosticsSeconds,
	double framePreparationSeconds,
	double frameWriteSeconds,
	double totalStepSeconds) {
	const fs::path path = diagnosticsDirectory() / "timing_detailed.csv";
	const bool shouldWriteHeader = !fs::exists(path);
	std::ofstream outFile = openDiagnosticsFile(path, std::ios::out | std::ios::app);
	const double measuredBackendSeconds = backendTiming.cellAssignmentSeconds
		+ backendTiming.neighborBuildSeconds
		+ backendTiming.wallActivationSeconds
		+ backendTiming.pressureSeconds
		+ backendTiming.wallUpdateSeconds
		+ backendTiming.wallResetSeconds
		+ backendTiming.rhsSeconds
		+ backendTiming.stateCopySeconds
		+ backendTiming.integrationSeconds
		+ backendTiming.stateApplySeconds
		+ backendTiming.cellClearSeconds;

	if (shouldWriteHeader) {
		outFile << "step,cell_assignment_seconds,neighbor_build_seconds,"
			<< "wall_activation_seconds,pressure_seconds,wall_update_seconds,"
			<< "wall_reset_seconds,rhs_seconds,state_copy_seconds,"
			<< "integration_seconds,state_apply_seconds,cell_clear_seconds,"
			<< "measured_backend_seconds,frame_diagnostics_seconds,"
			<< "frame_preparation_seconds,frame_write_seconds,total_step_seconds\n";
	}

	outFile << std::fixed << std::setprecision(12)
		<< step << ','
		<< backendTiming.cellAssignmentSeconds << ','
		<< backendTiming.neighborBuildSeconds << ','
		<< backendTiming.wallActivationSeconds << ','
		<< backendTiming.pressureSeconds << ','
		<< backendTiming.wallUpdateSeconds << ','
		<< backendTiming.wallResetSeconds << ','
		<< backendTiming.rhsSeconds << ','
		<< backendTiming.stateCopySeconds << ','
		<< backendTiming.integrationSeconds << ','
		<< backendTiming.stateApplySeconds << ','
		<< backendTiming.cellClearSeconds << ','
		<< measuredBackendSeconds << ','
		<< frameDiagnosticsSeconds << ','
		<< framePreparationSeconds << ','
		<< frameWriteSeconds << ','
		<< totalStepSeconds << '\n';
}
