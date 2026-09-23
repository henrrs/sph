#include "output_writer.h"

#include <stdexcept>

#include "legacy_vtk_writer.h"
#include "parameters.h"

#if defined(SPH_ENABLE_VTKHDF)
#include "vtkhdf_writer.h"
#endif

SimulationOutputWriter::SimulationOutputWriter() {
#if defined(SPH_ENABLE_VTKHDF)
	if constexpr (sim::outputFormat == sim::OutputFormat::VTKHDF || sim::outputFormat == sim::OutputFormat::Both) {
		vtkhdfWriter_ = std::make_unique<VtkHdfWriter>("simulation");
	}
#else
	if constexpr (sim::outputFormat == sim::OutputFormat::VTKHDF || sim::outputFormat == sim::OutputFormat::Both) {
		throw std::runtime_error("VTKHDF output was selected, but the binary was built without SPH_ENABLE_VTKHDF.");
	}
#endif
}

SimulationOutputWriter::~SimulationOutputWriter() = default;

void SimulationOutputWriter::writeFrame(const ParticleSet& particles,
	const ParticleSet& wallParticles,
	const std::string& name,
	int frameNumber,
	double physicalTime) {
	writeFrame(collectOutputFrameData(particles, wallParticles), name, frameNumber, physicalTime);
}

void SimulationOutputWriter::writeFrame(const OutputFrameData& frame,
	const std::string& name,
	int frameNumber,
	double physicalTime) {
	if constexpr (sim::outputFormat == sim::OutputFormat::LegacyVTK || sim::outputFormat == sim::OutputFormat::Both) {
		exportMeshLegacyVTK(frame, name, frameNumber);
	}

#if defined(SPH_ENABLE_VTKHDF)
	if constexpr (sim::outputFormat == sim::OutputFormat::VTKHDF || sim::outputFormat == sim::OutputFormat::Both) {
		vtkhdfWriter_->appendFrame(frame, physicalTime);
	}
#else
	(void)frame;
	(void)physicalTime;
#endif
}
