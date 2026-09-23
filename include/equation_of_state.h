#pragma once

#include <span>

double calcPressure(double densityOfParticleI);
double calcDensityFromPressure(double pressure);
void updatePressureField(std::span<const double> density, std::span<double> pressure);
