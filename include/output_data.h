#pragma once

#include <cstddef>
#include <vector>

#include "particle.h"

struct OutputFrameData {
	std::vector<double> points;
	std::vector<double> pressure;
	std::vector<double> velocity;
	std::vector<int> particleType;

	[[nodiscard]] std::size_t pointCount() const noexcept {
		return pressure.size();
	}
};

OutputFrameData collectOutputFrameData(const ParticleSet& particles, const ParticleSet& wallParticles);
