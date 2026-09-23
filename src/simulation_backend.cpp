#include "simulation_backend.h"

#include <memory>
#include <stdexcept>

#include "cpu_openmp_backend.h"

std::unique_ptr<SimulationBackend> createSimulationBackend(const SimulationConfig& config) {
	switch (config.backend) {
	case BackendKind::CpuOpenMp:
		return std::make_unique<CpuOpenMpBackend>(config);
	}

	throw std::runtime_error("Unsupported simulation backend.");
}
