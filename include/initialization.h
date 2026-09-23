#pragma once

#include <vector>

#include "cell.h"
#include "parameters.h"
#include "particle.h"

double calcInitialDensity(double y);

void startCells(std::vector<Cell>& cells);
void startParticles(ParticleSet& particles);
void startWalls(ParticleSet& wallParticles);
