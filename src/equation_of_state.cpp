#include "equation_of_state.h"

#include <algorithm>
#include <cmath>

#include "parameters.h"

namespace {

[[nodiscard]] inline double square(double value) noexcept {
	return value * value;
}

[[nodiscard]] inline double fourthPower(double value) noexcept {
	const double value2 = square(value);
	return value2 * value2;
}

[[nodiscard]] inline double seventhPower(double value) noexcept {
	return fourthPower(value) * square(value) * value;
}

[[nodiscard]] inline double equationOfStateB() noexcept {
	return (square(sim::speedOfSound) * sim::densityOfFluid) / sim::gamma;
}

} // namespace

double calcPressure(double densityOfParticleI) {
	static_assert(sim::gamma == 7.0, "calcPressure uses an optimized integer power for gamma = 7.");
	const double B = equationOfStateB();
	return B * (seventhPower(densityOfParticleI / sim::densityOfFluid) - 1.0);
}

double calcDensityFromPressure(double pressure) {
	const double B = equationOfStateB();
	const double densityRatio = std::max(1.0e-8, 1.0 + (pressure / B));
	return sim::densityOfFluid * std::pow(densityRatio, 1.0 / sim::gamma);
}

void updatePressureField(std::span<const double> density, std::span<double> pressure) {
	for (std::size_t i = 0; i < density.size(); i++) {
		pressure[i] = calcPressure(density[i]);
	}
}
