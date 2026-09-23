#include "particle.h"

void ParticleSet::resize(std::size_t count) {
	posX.resize(count);
	posY.resize(count);
	posZ.resize(count);
	velX.resize(count);
	velY.resize(count);
	velZ.resize(count);
	previousVelX.resize(count);
	previousVelY.resize(count);
	previousVelZ.resize(count);
	accX.resize(count);
	accY.resize(count);
	accZ.resize(count);
	xsphX.resize(count);
	xsphY.resize(count);
	xsphZ.resize(count);
	density.resize(count);
	densityDerivation.resize(count);
	pressure.resize(count);
	type.resize(count);
	tag.resize(count);
	currentCell.resize(count);
}

void ParticleSet::reserve(std::size_t count) {
	posX.reserve(count);
	posY.reserve(count);
	posZ.reserve(count);
	velX.reserve(count);
	velY.reserve(count);
	velZ.reserve(count);
	previousVelX.reserve(count);
	previousVelY.reserve(count);
	previousVelZ.reserve(count);
	accX.reserve(count);
	accY.reserve(count);
	accZ.reserve(count);
	xsphX.reserve(count);
	xsphY.reserve(count);
	xsphZ.reserve(count);
	density.reserve(count);
	densityDerivation.reserve(count);
	pressure.reserve(count);
	type.reserve(count);
	tag.reserve(count);
	currentCell.reserve(count);
}

void ParticleSet::clear() {
	posX.clear();
	posY.clear();
	posZ.clear();
	velX.clear();
	velY.clear();
	velZ.clear();
	previousVelX.clear();
	previousVelY.clear();
	previousVelZ.clear();
	accX.clear();
	accY.clear();
	accZ.clear();
	xsphX.clear();
	xsphY.clear();
	xsphZ.clear();
	density.clear();
	densityDerivation.clear();
	pressure.clear();
	type.clear();
	tag.clear();
	currentCell.clear();
}

std::size_t ParticleSet::size() const noexcept {
	return posX.size();
}

ParticleScalarsView makeParticleScalarsView(const ParticleSet& particles) noexcept {
	return {
		particles.posX,
		particles.posY,
		particles.posZ,
		particles.pressure,
		particles.velX,
		particles.velY,
		particles.velZ,
	};
}
