#pragma once

#include <cstddef>
#include <span>
#include <vector>

struct ParticleSet {
	std::vector<double> posX;
	std::vector<double> posY;
	std::vector<double> posZ;
	std::vector<double> velX;
	std::vector<double> velY;
	std::vector<double> velZ;
	std::vector<double> previousVelX;
	std::vector<double> previousVelY;
	std::vector<double> previousVelZ;
	std::vector<double> accX;
	std::vector<double> accY;
	std::vector<double> accZ;
	std::vector<double> xsphX;
	std::vector<double> xsphY;
	std::vector<double> xsphZ;
	std::vector<double> density;
	std::vector<double> densityDerivation;
	std::vector<double> pressure;
	std::vector<int> type;
	std::vector<int> tag;
	std::vector<int> currentCell;

	void resize(std::size_t count);
	void reserve(std::size_t count);
	void clear();
	[[nodiscard]] std::size_t size() const noexcept;
};

struct ParticleScalarsView {
	std::span<const double> posX;
	std::span<const double> posY;
	std::span<const double> posZ;
	std::span<const double> pressure;
	std::span<const double> velX;
	std::span<const double> velY;
	std::span<const double> velZ;
};

[[nodiscard]] ParticleScalarsView makeParticleScalarsView(const ParticleSet& particles) noexcept;
