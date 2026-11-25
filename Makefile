# Build the alphadual executable (C++17) with the vendored DAQP (C) and
# nanoflann (header-only).  Usage:  make            (optimised for this CPU)
#                                   make PORTABLE=1 (no -march=native)
# On Windows (MSYS2 / MinGW) the executable gets the .exe suffix.
EXE      := $(if $(filter Windows_NT,$(OS)),.exe,)
CXX      ?= g++
CC       ?= gcc
ARCH     := $(if $(PORTABLE),,-march=native)

# DAQP: same optimisation/math flags as its official Python package (which
# alpha_dualqp.py was validated against), but WITHOUT -DPROFILING, which
# would make DAQP read a clock inside every solve (internal timers).
DAQP_DIR := third_party/daqp
DAQP_FLAGS := -O3 $(ARCH) -fassociative-math -fno-signed-zeros -fno-trapping-math \
              -I$(DAQP_DIR)/include
CXXFLAGS := -O3 $(ARCH) -std=c++17 -Wall -Wextra -Wno-unused-parameter \
            -I$(DAQP_DIR)/include -Ithird_party/nanoflann

DAQP_SRC := $(wildcard $(DAQP_DIR)/src/*.c)
DAQP_OBJ := $(patsubst $(DAQP_DIR)/src/%.c,build/daqp_%.o,$(DAQP_SRC))

all: alphadual$(EXE)

alphadual$(EXE): build/alphadual.o $(DAQP_OBJ)
	$(CXX) -o $@ $^ -lm

build/alphadual.o: alphadual.cpp | build
	$(CXX) $(CXXFLAGS) -DALPHADUAL_BUILD_FLAGS='"$(strip $(CXXFLAGS))"' -c $< -o $@

build/daqp_%.o: $(DAQP_DIR)/src/%.c | build
	$(CC) $(DAQP_FLAGS) -c $< -o $@

build:
	mkdir -p build

# Development build with phase timers printed to stderr (not for benchmarks)
alphadual_profile$(EXE): alphadual.cpp $(DAQP_OBJ)
	$(CXX) $(CXXFLAGS) -DALPHADUAL_PROFILE -o $@ alphadual.cpp $(DAQP_OBJ) -lm

clean:
	rm -rf build alphadual alphadual_profile alphadual.exe alphadual_profile.exe

.PHONY: all clean
