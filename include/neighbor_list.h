#pragma once

#include <cstddef>
#include <vector>

struct NeighborListCSR {
	std::vector<int> fluidOffsets;
	std::vector<int> fluidIndices;
	std::vector<int> boundaryOffsets;
	std::vector<int> boundaryIndices;
	std::vector<int> wallFluidOffsets;
	std::vector<int> wallFluidIndices;
	std::vector<int> activeBoundaryIndices;

	void resizeOffsets(std::size_t particleCount) {
		fluidOffsets.assign(particleCount + 1, 0);
		boundaryOffsets.assign(particleCount + 1, 0);
	}

	void clear() {
		fluidOffsets.clear();
		fluidIndices.clear();
		boundaryOffsets.clear();
		boundaryIndices.clear();
		wallFluidOffsets.clear();
		wallFluidIndices.clear();
		activeBoundaryIndices.clear();
	}
};

struct NeighborBuildWorkspace {
	std::vector<int> fluidCounts;
	std::vector<int> boundaryCounts;
	std::vector<int> wallFluidCounts;
	std::vector<int> wallFluidWriteOffsets;
};
