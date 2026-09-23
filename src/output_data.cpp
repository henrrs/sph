#include "output_data.h"

#include "parameters.h"

namespace {

using namespace sim;

bool shouldExportWallParticle(const ParticleSet& wallParticles, std::size_t index) {
	if constexpr (!exportWallParticlesInVTK) {
		return false;
	}

	const double x = wallParticles.posX[index];
	const double y = wallParticles.posY[index];
	const double z = wallParticles.posZ[index];

	if constexpr (hideFrontWallParticlesInVTK) {
		if (z < 0.0) {
			return false;
		}
	}
	if constexpr (hideBackWallParticlesInVTK) {
		if (z > zDomain) {
			return false;
		}
	}
	if constexpr (hideLeftWallParticlesInVTK) {
		if (x < 0.0) {
			return false;
		}
	}
	if constexpr (hideRightWallParticlesInVTK) {
		if (x > xDomain) {
			return false;
		}
	}
	if constexpr (hideBottomWallParticlesInVTK) {
		if (y < 0.0) {
			return false;
		}
	}

	return true;
}

void appendParticle(OutputFrameData& frame, const ParticleSet& particles, std::size_t index) {
	frame.points.push_back(particles.posX[index]);
	frame.points.push_back(particles.posY[index]);
	frame.points.push_back(particles.posZ[index]);
	frame.pressure.push_back(particles.pressure[index]);
	frame.velocity.push_back(particles.velX[index]);
	frame.velocity.push_back(particles.velY[index]);
	frame.velocity.push_back(particles.velZ[index]);
	frame.particleType.push_back(particles.type[index]);
}

} // namespace

OutputFrameData collectOutputFrameData(const ParticleSet& particles, const ParticleSet& wallParticles) {
	OutputFrameData frame;
	frame.points.reserve((particles.size() + wallParticles.size()) * 3);
	frame.pressure.reserve(particles.size() + wallParticles.size());
	frame.velocity.reserve((particles.size() + wallParticles.size()) * 3);
	frame.particleType.reserve(particles.size() + wallParticles.size());

	for (std::size_t i = 0; i < particles.size(); i++) {
		appendParticle(frame, particles, i);
	}

	for (std::size_t i = 0; i < wallParticles.size(); i++) {
		if (shouldExportWallParticle(wallParticles, i)) {
			appendParticle(frame, wallParticles, i);
		}
	}

	return frame;
}
