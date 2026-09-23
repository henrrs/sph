#include "simulation_config.h"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

int parsePositiveInt(const std::string& value, const char* optionName) {
	const int parsed = std::stoi(value);
	if (parsed <= 0) {
		throw std::invalid_argument(std::string(optionName) + " must be positive.");
	}
	return parsed;
}

void printUsage(const char* executableName) {
	std::cout << "Usage: " << executableName << " [options]\n"
		<< "Options:\n"
		<< "  --steps N                    Override number of simulation steps.\n"
		<< "  --frames N                   Override number of output frames.\n"
		<< "  --backend cpu                Select CPU/OpenMP backend.\n"
		<< "  --disable-frame-diagnostics  Do not append output/diagnostics/simulation.csv.\n"
		<< "  --disable-timing             Do not append output/diagnostics/timing.csv.\n"
		<< "  --quiet                      Suppress per-step progress prints.\n"
		<< "  --help                       Show this message.\n";
}

} // namespace

SimulationConfig makeDefaultSimulationConfig() {
	return {};
}

SimulationConfig parseSimulationConfig(int argc, char** argv) {
	SimulationConfig config = makeDefaultSimulationConfig();

	for (int i = 1; i < argc; i++) {
		const std::string argument = argv[i];
		if (argument == "--help") {
			printUsage(argv[0]);
			std::exit(0);
		}
		else if (argument == "--steps") {
			if (i + 1 >= argc) {
				throw std::invalid_argument("--steps requires a value.");
			}
			config.numberOfSteps = parsePositiveInt(argv[++i], "--steps");
		}
		else if (argument == "--frames") {
			if (i + 1 >= argc) {
				throw std::invalid_argument("--frames requires a value.");
			}
			config.numFrames = parsePositiveInt(argv[++i], "--frames");
		}
		else if (argument == "--backend") {
			if (i + 1 >= argc) {
				throw std::invalid_argument("--backend requires a value.");
			}
			const std::string backend = argv[++i];
			if (backend == "cpu" || backend == "openmp") {
				config.backend = BackendKind::CpuOpenMp;
			}
			else {
				throw std::invalid_argument("Unsupported backend: " + backend);
			}
		}
		else if (argument == "--disable-frame-diagnostics") {
			config.enableFrameDiagnostics = false;
		}
		else if (argument == "--disable-timing") {
			config.enableTimingDiagnostics = false;
		}
		else if (argument == "--quiet") {
			config.printProgress = false;
		}
		else {
			throw std::invalid_argument("Unknown option: " + argument);
		}
	}

	if (config.numFrames > config.numberOfSteps) {
		config.numFrames = config.numberOfSteps;
	}

	return config;
}
