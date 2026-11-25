# Build gudhi_alpha (alpha filtration with GUDHI, any ambient dimension).
#   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release [-DGUDHI_INCLUDE_DIR=/path/to/gudhi/include]
#   cmake --build build --config Release
# The executable (gudhi_alpha, or gudhi_alpha.exe on Windows) is placed in THIS folder.
#
# Needs: GUDHI headers (shipped in third_party/gudhi, GUDHI 3.11.0, MIT license),
#        CGAL >= 5 (header-only), Eigen 3, Boost (headers), GMP and MPFR.
# Tell CMake where they are with -DCMAKE_PREFIX_PATH=... (e.g. a conda
# environment: "$CONDA_PREFIX/Library" on Windows, "$CONDA_PREFIX" elsewhere)
# and -DGUDHI_INCLUDE_DIR=... if GUDHI is not found automatically.
# Option -DGA_NATIVE=OFF: GCC/Clang without -march=native (portable binary).
cmake_minimum_required(VERSION 3.15)
project(gudhi_alpha LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
option(GA_NATIVE "GCC/Clang: optimise for this CPU (-march=native), as alphadual" ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(CGAL REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Boost REQUIRED)
find_path(GUDHI_INCLUDE_DIR gudhi/Alpha_complex.h
          HINTS ${CMAKE_CURRENT_SOURCE_DIR}/third_party/gudhi/include   # shipped (GUDHI 3.11.0)
                ENV GUDHI_INCLUDE_DIR PATH_SUFFIXES include)
if(NOT GUDHI_INCLUDE_DIR)
  message(FATAL_ERROR "GUDHI headers not found. Download a GUDHI release (https://github.com/GUDHI/gudhi-devel/releases) "
                      "and pass -DGUDHI_INCLUDE_DIR=<its include folder>")
endif()
message(STATUS "GUDHI headers: ${GUDHI_INCLUDE_DIR}")

add_executable(gudhi_alpha gudhi_alpha.cpp)
target_include_directories(gudhi_alpha PRIVATE ${GUDHI_INCLUDE_DIR})
target_link_libraries(gudhi_alpha PRIVATE CGAL::CGAL Eigen3::Eigen Boost::boost)
target_compile_definitions(gudhi_alpha PRIVATE CGAL_EIGEN3_ENABLED NOMINMAX)
if(WIN32)
  target_link_libraries(gudhi_alpha PRIVATE psapi)
endif()
if(MSVC)
  target_compile_options(gudhi_alpha PRIVATE /bigobj /permissive-)   # large CGAL templates
  target_compile_definitions(gudhi_alpha PRIVATE _CRT_SECURE_NO_WARNINGS _USE_MATH_DEFINES)
  set(GA_FLAGS_DESC "MSVC_Release_/O2")
else()
  target_compile_options(gudhi_alpha PRIVATE -O3 $<$<BOOL:${GA_NATIVE}>:-march=native>)
  if(GA_NATIVE)
    set(GA_FLAGS_DESC "-O3_-march=native")
  else()
    set(GA_FLAGS_DESC "-O3")
  endif()
endif()
target_compile_definitions(gudhi_alpha PRIVATE
  "GA_BUILD_FLAGS=\"cmake_${CMAKE_CXX_COMPILER_ID}_${GA_FLAGS_DESC}\"")

set_target_properties(gudhi_alpha PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
foreach(cfg Debug Release RelWithDebInfo MinSizeRel)
  string(TOUPPER ${cfg} CFG)
  set_target_properties(gudhi_alpha PROPERTIES RUNTIME_OUTPUT_DIRECTORY_${CFG} ${CMAKE_CURRENT_SOURCE_DIR})
endforeach()
