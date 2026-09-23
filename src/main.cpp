#include <exception>
#include <iostream>

#include "simulation.h"
#include "simulation_config.h"

int main(int argc, char** argv) {
	try {
		SimulationConfig config = parseSimulationConfig(argc, argv);
		SphSimulation simulation(config);
		simulation.initialize();
		simulation.run();
	}
	catch (const std::exception& exception) {
		std::cerr << "Error: " << exception.what() << '\n';
		return 1;
	}

	return 0;
}
