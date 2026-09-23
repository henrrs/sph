#pragma once

#include <cstddef>
#include <string>

#include "output_data.h"

class VtkHdfWriter {
public:
	explicit VtkHdfWriter(std::string fileStem);
	~VtkHdfWriter();

	VtkHdfWriter(const VtkHdfWriter&) = delete;
	VtkHdfWriter& operator=(const VtkHdfWriter&) = delete;

	void appendFrame(const OutputFrameData& frame, double physicalTime);

private:
	void initialize(std::size_t initialPointCount);
	void close() noexcept;

	std::string fileStem_;
	bool initialized_ = false;
	long long stepCount_ = 0;
	long long pointOffset_ = 0;
	long long partOffset_ = 0;
	long long vertexCellOffset_ = 0;
	long long vertexConnectivityOffset_ = 0;
	long long lineCellOffset_ = 0;
	long long lineConnectivityOffset_ = 0;
	long long polygonCellOffset_ = 0;
	long long polygonConnectivityOffset_ = 0;
	long long stripCellOffset_ = 0;
	long long stripConnectivityOffset_ = 0;

	long long fileHandle_ = -1;
};
