#pragma once

#include <memory>
#include <string>

#include "output_data.h"
#include "particle.h"

class VtkHdfWriter;

class SimulationOutputWriter {
public:
	SimulationOutputWriter();
	~SimulationOutputWriter();

	SimulationOutputWriter(const SimulationOutputWriter&) = delete;
	SimulationOutputWriter& operator=(const SimulationOutputWriter&) = delete;

	void writeFrame(const ParticleSet& particles,
		const ParticleSet& wallParticles,
		const std::string& name,
		int frameNumber,
		double physicalTime);
	void writeFrame(const OutputFrameData& frame,
		const std::string& name,
		int frameNumber,
		double physicalTime);

private:
#if defined(SPH_ENABLE_VTKHDF)
	std::unique_ptr<VtkHdfWriter> vtkhdfWriter_;
#endif
};
