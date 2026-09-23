#include "vtkhdf_writer.h"

#if !defined(SPH_ENABLE_VTKHDF)
#error "vtkhdf_writer.cpp must be compiled only when SPH_ENABLE_VTKHDF is defined."
#endif

#include <algorithm>
#include <filesystem>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include <hdf5.h>

#include "output_data.h"

namespace {

namespace fs = std::filesystem;

void checkStatus(herr_t status, const std::string& message) {
	if (status < 0) {
		throw std::runtime_error(message);
	}
}

hid_t checkedHandle(hid_t handle, const std::string& message) {
	if (handle < 0) {
		throw std::runtime_error(message);
	}
	return handle;
}

hid_t fileIdFromHandle(long long handle) {
	return static_cast<hid_t>(handle);
}

void closeIfValid(hid_t& handle, herr_t (*closeFn)(hid_t)) noexcept {
	if (handle >= 0) {
		closeFn(handle);
		handle = -1;
	}
}

hid_t createGroup(hid_t parent, const char* name) {
	return checkedHandle(H5Gcreate2(parent, name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
		std::string("Could not create HDF5 group: ") + name);
}

void createStringAttribute(hid_t object, const char* name, const std::string& value) {
	hid_t type = checkedHandle(H5Tcopy(H5T_C_S1), "Could not copy HDF5 string type.");
	checkStatus(H5Tset_size(type, value.size()), "Could not set HDF5 string attribute size.");
	checkStatus(H5Tset_strpad(type, H5T_STR_NULLPAD), "Could not set HDF5 string padding.");
	hid_t space = checkedHandle(H5Screate(H5S_SCALAR), "Could not create scalar dataspace for string attribute.");
	hid_t attribute = checkedHandle(H5Acreate2(object, name, type, space, H5P_DEFAULT, H5P_DEFAULT),
		std::string("Could not create HDF5 string attribute: ") + name);
	checkStatus(H5Awrite(attribute, type, value.c_str()), std::string("Could not write HDF5 string attribute: ") + name);
	closeIfValid(attribute, H5Aclose);
	closeIfValid(space, H5Sclose);
	closeIfValid(type, H5Tclose);
}

void writeLongLongScalarAttribute(hid_t object, const char* name, long long value) {
	hid_t space = checkedHandle(H5Screate(H5S_SCALAR), "Could not create scalar dataspace for integer attribute.");
	hid_t attribute = -1;
	const htri_t exists = H5Aexists(object, name);
	if (exists > 0) {
		attribute = checkedHandle(H5Aopen(object, name, H5P_DEFAULT), std::string("Could not open HDF5 attribute: ") + name);
	}
	else {
		attribute = checkedHandle(H5Acreate2(object, name, H5T_NATIVE_LLONG, space, H5P_DEFAULT, H5P_DEFAULT),
			std::string("Could not create HDF5 attribute: ") + name);
	}
	checkStatus(H5Awrite(attribute, H5T_NATIVE_LLONG, &value), std::string("Could not write HDF5 attribute: ") + name);
	closeIfValid(attribute, H5Aclose);
	closeIfValid(space, H5Sclose);
}

void createVersionAttribute(hid_t object) {
	const long long version[2] = { 1, 0 };
	const hsize_t dims[1] = { 2 };
	hid_t space = checkedHandle(H5Screate_simple(1, dims, nullptr), "Could not create Version dataspace.");
	hid_t attribute = checkedHandle(H5Acreate2(object, "Version", H5T_NATIVE_LLONG, space, H5P_DEFAULT, H5P_DEFAULT),
		"Could not create VTKHDF Version attribute.");
	checkStatus(H5Awrite(attribute, H5T_NATIVE_LLONG, version), "Could not write VTKHDF Version attribute.");
	closeIfValid(attribute, H5Aclose);
	closeIfValid(space, H5Sclose);
}

hid_t createResizableDataset(hid_t parent,
	const char* name,
	hid_t type,
	std::initializer_list<hsize_t> dims,
	std::initializer_list<hsize_t> maxDims,
	std::initializer_list<hsize_t> chunkDims) {
	std::vector<hsize_t> initial(dims);
	std::vector<hsize_t> maximum(maxDims);
	std::vector<hsize_t> chunk(chunkDims);

	hid_t space = checkedHandle(H5Screate_simple(static_cast<int>(initial.size()), initial.data(), maximum.data()),
		std::string("Could not create dataspace for dataset: ") + name);
	hid_t creationProperties = checkedHandle(H5Pcreate(H5P_DATASET_CREATE),
		std::string("Could not create dataset properties for: ") + name);
	checkStatus(H5Pset_chunk(creationProperties, static_cast<int>(chunk.size()), chunk.data()),
		std::string("Could not set chunking for dataset: ") + name);
	hid_t dataset = checkedHandle(H5Dcreate2(parent, name, type, space, H5P_DEFAULT, creationProperties, H5P_DEFAULT),
		std::string("Could not create HDF5 dataset: ") + name);
	closeIfValid(creationProperties, H5Pclose);
	closeIfValid(space, H5Sclose);
	return dataset;
}

void createOneDimensionalDataset(hid_t parent, const char* name, hid_t type, hsize_t chunkSize) {
	hid_t dataset = createResizableDataset(parent, name, type, { 0 }, { H5S_UNLIMITED }, { std::max<hsize_t>(1, chunkSize) });
	closeIfValid(dataset, H5Dclose);
}

void createTwoDimensionalDataset(hid_t parent, const char* name, hid_t type, hsize_t chunkRows, hsize_t componentCount) {
	hid_t dataset = createResizableDataset(parent,
		name,
		type,
		{ 0, componentCount },
		{ H5S_UNLIMITED, componentCount },
		{ std::max<hsize_t>(1, chunkRows), componentCount });
	closeIfValid(dataset, H5Dclose);
}

void append1D(hid_t file, const char* path, hid_t type, const void* values, hsize_t valueCount) {
	if (valueCount == 0) {
		return;
	}

	hid_t dataset = checkedHandle(H5Dopen2(file, path, H5P_DEFAULT), std::string("Could not open dataset: ") + path);
	hid_t oldSpace = checkedHandle(H5Dget_space(dataset), std::string("Could not get dataspace for dataset: ") + path);
	hsize_t oldDims[1] = { 0 };
	checkStatus(H5Sget_simple_extent_dims(oldSpace, oldDims, nullptr), std::string("Could not read dataset dimensions: ") + path);
	closeIfValid(oldSpace, H5Sclose);

	const hsize_t newDims[1] = { oldDims[0] + valueCount };
	checkStatus(H5Dset_extent(dataset, newDims), std::string("Could not resize dataset: ") + path);
	hid_t fileSpace = checkedHandle(H5Dget_space(dataset), std::string("Could not get resized dataspace: ") + path);
	const hsize_t start[1] = { oldDims[0] };
	const hsize_t count[1] = { valueCount };
	checkStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, start, nullptr, count, nullptr),
		std::string("Could not select hyperslab for dataset: ") + path);
	hid_t memorySpace = checkedHandle(H5Screate_simple(1, count, nullptr), std::string("Could not create memory dataspace: ") + path);
	checkStatus(H5Dwrite(dataset, type, memorySpace, fileSpace, H5P_DEFAULT, values),
		std::string("Could not append data to dataset: ") + path);
	closeIfValid(memorySpace, H5Sclose);
	closeIfValid(fileSpace, H5Sclose);
	closeIfValid(dataset, H5Dclose);
}

void append2D(hid_t file, const char* path, hid_t type, const void* values, hsize_t rowCount, hsize_t componentCount) {
	if (rowCount == 0) {
		return;
	}

	hid_t dataset = checkedHandle(H5Dopen2(file, path, H5P_DEFAULT), std::string("Could not open dataset: ") + path);
	hid_t oldSpace = checkedHandle(H5Dget_space(dataset), std::string("Could not get dataspace for dataset: ") + path);
	hsize_t oldDims[2] = { 0, componentCount };
	checkStatus(H5Sget_simple_extent_dims(oldSpace, oldDims, nullptr), std::string("Could not read dataset dimensions: ") + path);
	closeIfValid(oldSpace, H5Sclose);

	const hsize_t newDims[2] = { oldDims[0] + rowCount, componentCount };
	checkStatus(H5Dset_extent(dataset, newDims), std::string("Could not resize dataset: ") + path);
	hid_t fileSpace = checkedHandle(H5Dget_space(dataset), std::string("Could not get resized dataspace: ") + path);
	const hsize_t start[2] = { oldDims[0], 0 };
	const hsize_t count[2] = { rowCount, componentCount };
	checkStatus(H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, start, nullptr, count, nullptr),
		std::string("Could not select hyperslab for dataset: ") + path);
	hid_t memorySpace = checkedHandle(H5Screate_simple(2, count, nullptr), std::string("Could not create memory dataspace: ") + path);
	checkStatus(H5Dwrite(dataset, type, memorySpace, fileSpace, H5P_DEFAULT, values),
		std::string("Could not append data to dataset: ") + path);
	closeIfValid(memorySpace, H5Sclose);
	closeIfValid(fileSpace, H5Sclose);
	closeIfValid(dataset, H5Dclose);
}

void createTopologyGroup(hid_t vtkhdfGroup, const char* name, hsize_t pointChunkSize) {
	hid_t group = createGroup(vtkhdfGroup, name);
	createOneDimensionalDataset(group, "Connectivity", H5T_NATIVE_LLONG, pointChunkSize);
	createOneDimensionalDataset(group, "Offsets", H5T_NATIVE_LLONG, pointChunkSize + 1);
	createOneDimensionalDataset(group, "NumberOfCells", H5T_NATIVE_LLONG, 128);
	createOneDimensionalDataset(group, "NumberOfConnectivityIds", H5T_NATIVE_LLONG, 128);
	closeIfValid(group, H5Gclose);
}

void createPointDataDataset(hid_t pointDataGroup,
	const char* name,
	hid_t type,
	hsize_t pointChunkSize,
	hsize_t componentCount,
	const char* attributeRole) {
	hid_t dataset = -1;
	if (componentCount == 1) {
		dataset = createResizableDataset(pointDataGroup, name, type, { 0 }, { H5S_UNLIMITED }, { std::max<hsize_t>(1, pointChunkSize) });
	}
	else {
		dataset = createResizableDataset(pointDataGroup,
			name,
			type,
			{ 0, componentCount },
			{ H5S_UNLIMITED, componentCount },
			{ std::max<hsize_t>(1, pointChunkSize), componentCount });
	}
	createStringAttribute(dataset, "Attribute", attributeRole);
	closeIfValid(dataset, H5Dclose);
}

void appendTopologyFrame(hid_t file,
	const char* groupName,
	long long numberOfCells,
	long long numberOfConnectivityIds,
	const std::vector<long long>& connectivity,
	const std::vector<long long>& offsets) {
	const std::string basePath = std::string("/VTKHDF/") + groupName;
	append1D(file, (basePath + "/NumberOfCells").c_str(), H5T_NATIVE_LLONG, &numberOfCells, 1);
	append1D(file, (basePath + "/NumberOfConnectivityIds").c_str(), H5T_NATIVE_LLONG, &numberOfConnectivityIds, 1);
	append1D(file, (basePath + "/Connectivity").c_str(), H5T_NATIVE_LLONG, connectivity.data(), connectivity.size());
	append1D(file, (basePath + "/Offsets").c_str(), H5T_NATIVE_LLONG, offsets.data(), offsets.size());
}

std::vector<long long> makeVertexConnectivity(std::size_t pointCount) {
	std::vector<long long> connectivity(pointCount);
	std::iota(connectivity.begin(), connectivity.end(), 0LL);
	return connectivity;
}

std::vector<long long> makeVertexOffsets(std::size_t pointCount) {
	std::vector<long long> offsets(pointCount + 1);
	std::iota(offsets.begin(), offsets.end(), 0LL);
	return offsets;
}

} // namespace

VtkHdfWriter::VtkHdfWriter(std::string fileStem)
	: fileStem_(std::move(fileStem)) {}

VtkHdfWriter::~VtkHdfWriter() {
	close();
}

void VtkHdfWriter::initialize(std::size_t initialPointCount) {
	fs::create_directories("output/vtkhdf");
	const fs::path filePath = fs::path("output/vtkhdf") / (fileStem_ + ".vtkhdf");
	const hsize_t pointChunkSize = std::clamp<hsize_t>(static_cast<hsize_t>(initialPointCount), 1, 65536);

	hid_t file = checkedHandle(H5Fcreate(filePath.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT),
		"Could not create VTKHDF file: " + filePath.string());
	fileHandle_ = static_cast<long long>(file);

	hid_t vtkhdfGroup = createGroup(file, "VTKHDF");
	createStringAttribute(vtkhdfGroup, "Type", "PolyData");
	createVersionAttribute(vtkhdfGroup);

	createOneDimensionalDataset(vtkhdfGroup, "NumberOfPoints", H5T_NATIVE_LLONG, 128);
	createTwoDimensionalDataset(vtkhdfGroup, "Points", H5T_NATIVE_DOUBLE, pointChunkSize, 3);

	hid_t fieldDataGroup = createGroup(vtkhdfGroup, "FieldData");
	closeIfValid(fieldDataGroup, H5Gclose);
	hid_t cellDataGroup = createGroup(vtkhdfGroup, "CellData");
	closeIfValid(cellDataGroup, H5Gclose);

	hid_t pointDataGroup = createGroup(vtkhdfGroup, "PointData");
	createStringAttribute(pointDataGroup, "Scalars", "pressure");
	createStringAttribute(pointDataGroup, "Vectors", "velocity");
	createPointDataDataset(pointDataGroup, "pressure", H5T_NATIVE_DOUBLE, pointChunkSize, 1, "Scalars");
	createPointDataDataset(pointDataGroup, "velocity", H5T_NATIVE_DOUBLE, pointChunkSize, 3, "Vectors");
	createPointDataDataset(pointDataGroup, "particle_type", H5T_NATIVE_INT, pointChunkSize, 1, "Scalars");
	closeIfValid(pointDataGroup, H5Gclose);

	createTopologyGroup(vtkhdfGroup, "Vertices", pointChunkSize);
	createTopologyGroup(vtkhdfGroup, "Lines", 128);
	createTopologyGroup(vtkhdfGroup, "Polygons", 128);
	createTopologyGroup(vtkhdfGroup, "Strips", 128);

	hid_t stepsGroup = createGroup(vtkhdfGroup, "Steps");
	writeLongLongScalarAttribute(stepsGroup, "NSteps", 0);
	createOneDimensionalDataset(stepsGroup, "Values", H5T_NATIVE_DOUBLE, 128);
	createOneDimensionalDataset(stepsGroup, "PartOffsets", H5T_NATIVE_LLONG, 128);
	createOneDimensionalDataset(stepsGroup, "NumberOfParts", H5T_NATIVE_LLONG, 128);
	createOneDimensionalDataset(stepsGroup, "PointOffsets", H5T_NATIVE_LLONG, 128);
	createTwoDimensionalDataset(stepsGroup, "CellOffsets", H5T_NATIVE_LLONG, 128, 4);
	createTwoDimensionalDataset(stepsGroup, "ConnectivityIdOffsets", H5T_NATIVE_LLONG, 128, 4);

	hid_t pointDataOffsetsGroup = createGroup(stepsGroup, "PointDataOffsets");
	createOneDimensionalDataset(pointDataOffsetsGroup, "pressure", H5T_NATIVE_LLONG, 128);
	createOneDimensionalDataset(pointDataOffsetsGroup, "velocity", H5T_NATIVE_LLONG, 128);
	createOneDimensionalDataset(pointDataOffsetsGroup, "particle_type", H5T_NATIVE_LLONG, 128);
	closeIfValid(pointDataOffsetsGroup, H5Gclose);
	closeIfValid(stepsGroup, H5Gclose);
	closeIfValid(vtkhdfGroup, H5Gclose);

	initialized_ = true;
}

void VtkHdfWriter::appendFrame(const OutputFrameData& frame, double physicalTime) {
	if (!initialized_) {
		initialize(frame.pointCount());
	}

	hid_t file = fileIdFromHandle(fileHandle_);
	const long long pointCount = static_cast<long long>(frame.pointCount());
	const std::vector<long long> vertexConnectivity = makeVertexConnectivity(frame.pointCount());
	const std::vector<long long> vertexOffsets = makeVertexOffsets(frame.pointCount());
	const std::vector<long long> emptyOffsets = { 0 };
	const std::vector<long long> emptyConnectivity;

	append1D(file, "/VTKHDF/NumberOfPoints", H5T_NATIVE_LLONG, &pointCount, 1);
	append2D(file, "/VTKHDF/Points", H5T_NATIVE_DOUBLE, frame.points.data(), frame.pointCount(), 3);
	append1D(file, "/VTKHDF/PointData/pressure", H5T_NATIVE_DOUBLE, frame.pressure.data(), frame.pressure.size());
	append2D(file, "/VTKHDF/PointData/velocity", H5T_NATIVE_DOUBLE, frame.velocity.data(), frame.pointCount(), 3);
	append1D(file, "/VTKHDF/PointData/particle_type", H5T_NATIVE_INT, frame.particleType.data(), frame.particleType.size());

	appendTopologyFrame(file, "Vertices", pointCount, pointCount, vertexConnectivity, vertexOffsets);
	appendTopologyFrame(file, "Lines", 0, 0, emptyConnectivity, emptyOffsets);
	appendTopologyFrame(file, "Polygons", 0, 0, emptyConnectivity, emptyOffsets);
	appendTopologyFrame(file, "Strips", 0, 0, emptyConnectivity, emptyOffsets);

	const long long cellOffsets[4] = {
		vertexCellOffset_,
		lineCellOffset_,
		polygonCellOffset_,
		stripCellOffset_,
	};
	const long long connectivityIdOffsets[4] = {
		vertexConnectivityOffset_,
		lineConnectivityOffset_,
		polygonConnectivityOffset_,
		stripConnectivityOffset_,
	};
	const long long numberOfParts = 1;

	append1D(file, "/VTKHDF/Steps/Values", H5T_NATIVE_DOUBLE, &physicalTime, 1);
	append1D(file, "/VTKHDF/Steps/PartOffsets", H5T_NATIVE_LLONG, &partOffset_, 1);
	append1D(file, "/VTKHDF/Steps/NumberOfParts", H5T_NATIVE_LLONG, &numberOfParts, 1);
	append1D(file, "/VTKHDF/Steps/PointOffsets", H5T_NATIVE_LLONG, &pointOffset_, 1);
	append2D(file, "/VTKHDF/Steps/CellOffsets", H5T_NATIVE_LLONG, cellOffsets, 1, 4);
	append2D(file, "/VTKHDF/Steps/ConnectivityIdOffsets", H5T_NATIVE_LLONG, connectivityIdOffsets, 1, 4);
	append1D(file, "/VTKHDF/Steps/PointDataOffsets/pressure", H5T_NATIVE_LLONG, &pointOffset_, 1);
	append1D(file, "/VTKHDF/Steps/PointDataOffsets/velocity", H5T_NATIVE_LLONG, &pointOffset_, 1);
	append1D(file, "/VTKHDF/Steps/PointDataOffsets/particle_type", H5T_NATIVE_LLONG, &pointOffset_, 1);

	pointOffset_ += pointCount;
	partOffset_ += 1;
	vertexCellOffset_ += pointCount;
	vertexConnectivityOffset_ += pointCount;
	stepCount_ += 1;

	hid_t stepsGroup = checkedHandle(H5Gopen2(file, "/VTKHDF/Steps", H5P_DEFAULT), "Could not open /VTKHDF/Steps.");
	writeLongLongScalarAttribute(stepsGroup, "NSteps", stepCount_);
	closeIfValid(stepsGroup, H5Gclose);
	checkStatus(H5Fflush(file, H5F_SCOPE_GLOBAL), "Could not flush VTKHDF file.");
}

void VtkHdfWriter::close() noexcept {
	if (fileHandle_ >= 0) {
		hid_t file = fileIdFromHandle(fileHandle_);
		H5Fflush(file, H5F_SCOPE_GLOBAL);
		H5Fclose(file);
		fileHandle_ = -1;
	}
}
