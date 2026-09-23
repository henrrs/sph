#pragma once

#include <string>
#include <vector>

#include "output_data.h"
#include "particle.h"

void createFile(const ParticleSet& particles, std::string name, int frameNumber);
void createFilePositionMomentum(const ParticleSet& particles, std::string name, int frameNumber);
void createPositionFile(const ParticleSet& particles, std::string name, int frameNumber);
void createMomentumFile(const ParticleSet& particles, std::string name, int frameNumber);
void exportMeshLegacyVTK(const OutputFrameData& frame, std::string name, int frameNumber);
