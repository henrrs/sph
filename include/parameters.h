#pragma once

#include <cstddef>

namespace sim {

constexpr int ceilToInt(double value) {
	const int truncated = static_cast<int>(value);
	return (value > static_cast<double>(truncated)) ? truncated + 1 : truncated;
}

enum class OutputFormat {
	LegacyVTK,
	VTKHDF,
	Both,
};

#if defined(SPH_OUTPUT_FORMAT_LEGACY)
inline constexpr OutputFormat outputFormat = OutputFormat::LegacyVTK;
#elif defined(SPH_OUTPUT_FORMAT_BOTH)
inline constexpr OutputFormat outputFormat = OutputFormat::Both;
#else
inline constexpr OutputFormat outputFormat = OutputFormat::VTKHDF;
#endif

inline constexpr int numFrames = 500;
inline constexpr int numberOfSteps = 30000;
inline constexpr bool exportWallParticlesInVTK = true;
inline constexpr bool hideFrontWallParticlesInVTK = true;
inline constexpr bool hideBackWallParticlesInVTK = true;
inline constexpr bool hideLeftWallParticlesInVTK = false;
inline constexpr bool hideRightWallParticlesInVTK = false;
inline constexpr bool hideBottomWallParticlesInVTK = false;

inline constexpr double alpha = -0.40;
inline constexpr double beta = 0.70;
inline constexpr double gamma = 7.0;

inline constexpr double xsphConstant = -0.001;

inline constexpr double densityOfFluid = 1000.0;
inline constexpr double gravity = 9.81;
inline constexpr double pi = 3.14159265359;
inline constexpr double dt = 0.00005;
inline constexpr bool useVelocityVerletIntegrator = true;

inline constexpr int lineParticles = 20;
inline constexpr int columnParticles = 40;
inline constexpr int depthParticles = 20;
inline constexpr int numFluidParticles = lineParticles * columnParticles * depthParticles;

inline constexpr double widthCell = 20.0 / 100.0;
inline constexpr double heightCell = 40.0 / 100.0;
inline constexpr double depthCell = 20.0 / 100.0;

inline constexpr int xTimesDomain = 5;
inline constexpr int yTimesDomain = 1;

inline constexpr double xDomain = xTimesDomain * widthCell;
inline constexpr double yDomain = yTimesDomain * heightCell;
inline constexpr double zDomain = depthCell;

inline constexpr double space = 0.01;
inline constexpr double h = 1.2 * space;
inline constexpr int supportRadiusFactor = 3;
inline constexpr double kh = supportRadiusFactor * h;

inline constexpr int wallLayerCount = 6;
inline constexpr double wallThickness = wallLayerCount * space;
inline constexpr int minimumWallLayerCount = ceilToInt(kh / space);
static_assert(wallLayerCount >= minimumWallLayerCount,
	"The wall must be at least as thick as the kernel support radius.");

inline constexpr double domainMinX = -wallThickness;
inline constexpr double domainMaxX = xDomain + wallThickness;
inline constexpr double domainMinY = -wallThickness;
inline constexpr double domainMaxY = yDomain + wallThickness;
inline constexpr double domainMinZ = -wallThickness;
inline constexpr double domainMaxZ = zDomain + wallThickness;

inline constexpr int numberOfCellsInX = ceilToInt((domainMaxX - domainMinX) / kh);
inline constexpr int numberOfCellsInY = ceilToInt((domainMaxY - domainMinY) / kh);
inline constexpr int numberOfCellsInZ = ceilToInt((domainMaxZ - domainMinZ) / kh);
inline constexpr int totalNumCells = numberOfCellsInX * numberOfCellsInY * numberOfCellsInZ;

inline constexpr double speedOfSound = 28.0142820707074;
inline constexpr double particleMass = densityOfFluid * space * space * space;
inline constexpr double referenceParticleVolume = particleMass / densityOfFluid;

inline constexpr int wallParticleType = 1;
inline constexpr int fluidParticleType = 2;
inline constexpr double fluidHeight = yDomain;

} // namespace sim
