#pragma once

#include "particle.h"
#include "simulation_backend.h"

void writeHydrostaticInitializationDiagnostics(const ParticleSet& particles);
void appendSimulationDiagnostics(const ParticleSet& particles,
	int step,
	double physicalTime,
	int activeWallParticleCount,
	int fluidNeighborPairCount,
	int boundaryNeighborPairCount);
void appendTimingDiagnostics(int step,
	double computeSeconds,
	double frameSeconds,
	double totalSeconds);
void appendDetailedTimingDiagnostics(int step,
	const BackendStepTiming& backendTiming,
	double frameDiagnosticsSeconds,
	double framePreparationSeconds,
	double frameWriteSeconds,
	double totalStepSeconds);
