#pragma once

#include <array>
#include <vector>

struct Cell {
	double x0 = 0.0;
	double x1 = 0.0;
	double y0 = 0.0;
	double y1 = 0.0;
	double z0 = 0.0;
	double z1 = 0.0;

	int tag = 0;
	int line = 0;
	int column = 0;
	int depth = 0;

	std::vector<int> particles;
	std::vector<int> boundaryParticles;
	std::array<int, 27> neighborTags{};
	int neighborCount = 0;

	void clearParticles();
};
