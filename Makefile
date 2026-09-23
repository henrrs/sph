SHELL := /bin/bash

PROFILE ?= release
OUTPUT_FORMAT ?= vtkhdf

CXX ?= g++
OPENMP_FLAGS := -fopenmp
CPPFLAGS ?= -Iinclude
BASE_CXXFLAGS := -std=c++20 -Wall -Wextra -Wpedantic $(OPENMP_FLAGS)
LDFLAGS ?=
LDLIBS ?= $(OPENMP_FLAGS)

ifeq ($(PROFILE),default)
OPTFLAGS :=
else ifeq ($(PROFILE),debug)
OPTFLAGS := -O0 -g
else ifeq ($(PROFILE),release)
OPTFLAGS := -O2
else ifeq ($(PROFILE),o2)
OPTFLAGS := -O2
else ifeq ($(PROFILE),o3)
OPTFLAGS := -O3
else ifeq ($(PROFILE),ofast)
OPTFLAGS := -Ofast
else
$(error Unknown PROFILE '$(PROFILE)'. Use default, debug, release, o2, o3, or ofast)
endif

ifeq ($(OUTPUT_FORMAT),legacy)
OUTPUT_FORMAT_FLAGS := -DSPH_OUTPUT_FORMAT_LEGACY
REQUIRES_HDF5 := 0
else ifeq ($(OUTPUT_FORMAT),vtkhdf)
OUTPUT_FORMAT_FLAGS := -DSPH_OUTPUT_FORMAT_VTKHDF
REQUIRES_HDF5 := 1
else ifeq ($(OUTPUT_FORMAT),both)
OUTPUT_FORMAT_FLAGS := -DSPH_OUTPUT_FORMAT_BOTH
REQUIRES_HDF5 := 1
else
$(error Unknown OUTPUT_FORMAT '$(OUTPUT_FORMAT)'. Use legacy, vtkhdf, or both)
endif

ifeq ($(REQUIRES_HDF5),1)
HDF5_CFLAGS := $(shell pkg-config --cflags hdf5 2>/dev/null)
HDF5_LIBS := $(shell pkg-config --libs hdf5 2>/dev/null)
ifeq ($(strip $(HDF5_LIBS)),)
$(error OUTPUT_FORMAT=$(OUTPUT_FORMAT) requires HDF5. Install it with 'sudo apt install libhdf5-dev' or build with OUTPUT_FORMAT=legacy)
endif
CPPFLAGS += -DSPH_ENABLE_VTKHDF $(HDF5_CFLAGS)
LDLIBS += $(HDF5_LIBS)
endif

CPPFLAGS += $(OUTPUT_FORMAT_FLAGS)
CXXFLAGS ?= $(BASE_CXXFLAGS) $(OPTFLAGS)

BUILD_ROOT := build/$(PROFILE)
BUILD_DIR := $(BUILD_ROOT)/$(OUTPUT_FORMAT)
TARGET := $(BUILD_DIR)/sph
TARGET_ALIAS := $(BUILD_ROOT)/sph
TEST_BUILD_DIR := $(BUILD_ROOT)/tests
PHYSICS_TEST_TARGET := $(TEST_BUILD_DIR)/physics_consistency

SOURCES := \
	src/main.cpp \
	src/cell.cpp \
	src/cpu_openmp_backend.cpp \
	src/diagnostics.cpp \
	src/equation_of_state.cpp \
	src/initialization.cpp \
	src/integrators.cpp \
	src/legacy_vtk_writer.cpp \
	src/neighborhood.cpp \
	src/output_data.cpp \
	src/output_writer.cpp \
	src/particle.cpp \
	src/simulation.cpp \
	src/simulation_backend.cpp \
	src/simulation_config.cpp \
	src/sph_kernels.cpp \
	src/sph_rhs.cpp \
	src/wall_boundary.cpp

ifeq ($(REQUIRES_HDF5),1)
SOURCES += src/vtkhdf_writer.cpp
endif

OBJECTS := $(patsubst src/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))
DEPFILES := $(OBJECTS:.o=.d)
PHYSICS_TEST_SOURCES := \
	tests/physics_consistency.cpp \
	src/cell.cpp \
	src/equation_of_state.cpp \
	src/initialization.cpp \
	src/integrators.cpp \
	src/neighborhood.cpp \
	src/particle.cpp \
	src/sph_kernels.cpp \
	src/sph_rhs.cpp \
	src/wall_boundary.cpp

.PHONY: all default debug release o2 o3 ofast bench run test test-physics clean clean-output print-config FORCE

all: $(TARGET_ALIAS)

$(TARGET): $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)

$(TARGET_ALIAS): $(TARGET) FORCE
	ln -sfn $(OUTPUT_FORMAT)/sph $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

output:
	mkdir -p output/diagnostics output/mat output/vtk output/vtkhdf

$(BUILD_DIR)/%.o: src/%.cpp | $(BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

run: $(TARGET) | output
	./$(TARGET)

$(TEST_BUILD_DIR):
	mkdir -p $(TEST_BUILD_DIR)

$(PHYSICS_TEST_TARGET): $(PHYSICS_TEST_SOURCES) $(wildcard include/*.h) | $(TEST_BUILD_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(PHYSICS_TEST_SOURCES) -o $@ $(LDLIBS)

test: test-physics

test-physics: $(PHYSICS_TEST_TARGET)
	OMP_NUM_THREADS="$${OMP_NUM_THREADS:-$$(nproc)}" ./$(PHYSICS_TEST_TARGET)

default:
	$(MAKE) PROFILE=default all

debug:
	$(MAKE) PROFILE=debug all

release:
	$(MAKE) PROFILE=release all

o2:
	$(MAKE) PROFILE=o2 all

o3:
	$(MAKE) PROFILE=o3 all

ofast:
	$(MAKE) PROFILE=ofast all

bench:
	mkdir -p docs
	if [[ ! -f docs/benchmark.md ]]; then \
		echo "# Benchmark" > docs/benchmark.md; \
		echo "" >> docs/benchmark.md; \
	fi
	echo "" >> docs/benchmark.md
	echo "## Benchmark completo" >> docs/benchmark.md
	echo "" >> docs/benchmark.md
	echo "Data: $$(date '+%Y-%m-%d %H:%M:%S %Z')" >> docs/benchmark.md
	echo "" >> docs/benchmark.md
	echo "Este benchmark executa a simulacao completa em quatro perfis de compilacao, limpando \`output/\` entre as rodadas." >> docs/benchmark.md
	echo "" >> docs/benchmark.md
	echo "| Perfil | Flags | Tempo (s) | Speedup vs default | Speedup vs O2 |" >> docs/benchmark.md
	echo "| --- | --- | ---: | ---: | ---: |" >> docs/benchmark.md
	default_time=""; o2_time=""; \
	for profile in default o2 o3 ofast; do \
		case "$$profile" in \
			default) flags="(compiler default, sem -O explicito)" ;; \
			o2) flags="-O2" ;; \
			o3) flags="-O3" ;; \
			ofast) flags="-Ofast" ;; \
		esac; \
		echo "==> Benchmark $$profile ($$flags)"; \
		$(MAKE) clean-output >/dev/null; \
		$(MAKE) PROFILE=$$profile all >/dev/null; \
		/usr/bin/time -f "%e" -o .bench_time ./build/$$profile/sph >/dev/null 2>&1; \
		time_value="$$(cat .bench_time)"; \
		if [[ "$$profile" == "default" ]]; then default_time="$$time_value"; fi; \
		if [[ "$$profile" == "o2" ]]; then o2_time="$$time_value"; fi; \
		speedup_default="-"; \
		speedup_o2="-"; \
		if [[ -n "$$default_time" ]]; then speedup_default="$$(awk -v base="$$default_time" -v current="$$time_value" 'BEGIN { printf("%.3f", base / current) }')x"; fi; \
		if [[ -n "$$o2_time" ]]; then speedup_o2="$$(awk -v base="$$o2_time" -v current="$$time_value" 'BEGIN { printf("%.3f", base / current) }')x"; fi; \
		echo "| $$profile | $$flags | $$time_value | $$speedup_default | $$speedup_o2 |" >> docs/benchmark.md; \
		$(MAKE) clean-output >/dev/null; \
	done; \
	rm -f .bench_time

print-config:
	@echo PROFILE=$(PROFILE)
	@echo OUTPUT_FORMAT=$(OUTPUT_FORMAT)
	@echo REQUIRES_HDF5=$(REQUIRES_HDF5)
	@echo BUILD_DIR=$(BUILD_DIR)
	@echo TARGET=$(TARGET)
	@echo TARGET_ALIAS=$(TARGET_ALIAS)
	@echo CXXFLAGS=$(CXXFLAGS)
	@echo CPPFLAGS=$(CPPFLAGS)

clean:
	rm -rf build

clean-output:
	rm -f output/vtk/*.vtk output/vtkhdf/*.vtkhdf output/mat/*.mat output/diagnostics/*.csv

-include $(DEPFILES)
