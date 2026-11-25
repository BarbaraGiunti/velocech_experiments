# Build alphadual with CMake (alternative to the Makefile), with any compiler
# CMake finds: GCC, Clang, MinGW, or Visual C++ (MSVC). From this folder:
#   cmake -S . -B build-cmake -DCMAKE_BUILD_TYPE=Release
#   cmake --build build-cmake --config Release
# The executable (alphadual, or alphadual.exe on Windows) is placed in THIS
# folder, next to the Python and bash scripts that run it.
# Options: -DALPHADUAL_NATIVE=OFF  (GCC/Clang: no -march=native, portable binary)
cmake_minimum_required(VERSION 3.15)
project(alphadual LANGUAGES C CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
option(ALPHADUAL_NATIVE "GCC/Clang: optimise for this CPU (-march=native)" ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB DAQP_SRC ${CMAKE_CURRENT_SOURCE_DIR}/third_party/daqp/src/*.c)
add_executable(alphadual alphadual.cpp ${DAQP_SRC})
target_include_directories(alphadual PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/daqp/include
  ${CMAKE_CURRENT_SOURCE_DIR}/third_party/nanoflann)
find_package(Threads REQUIRED)
target_link_libraries(alphadual PRIVATE Threads::Threads)
target_compile_definitions(alphadual PRIVATE NOMINMAX)

if(MSVC)
  # Release already uses /O2. No DAQP PROFILING define: no timers inside DAQP.
  target_compile_definitions(alphadual PRIVATE _CRT_SECURE_NO_WARNINGS)
  set(ALPHADUAL_FLAGS_DESC "MSVC_Release_/O2")
else()
  target_link_libraries(alphadual PRIVATE m)
  target_compile_options(alphadual PRIVATE -O3
    $<$<BOOL:${ALPHADUAL_NATIVE}>:-march=native>)
  # DAQP: same math flags as its official Python package (as in the Makefile)
  set_source_files_properties(${DAQP_SRC} PROPERTIES COMPILE_OPTIONS
    "-fassociative-math;-fno-signed-zeros;-fno-trapping-math")
  if(ALPHADUAL_NATIVE)
    set(ALPHADUAL_FLAGS_DESC "-O3_-march=native")
  else()
    set(ALPHADUAL_FLAGS_DESC "-O3")
  endif()
endif()
# Recorded by "alphadual --version" (and so in the Excel report).
target_compile_definitions(alphadual PRIVATE
  "ALPHADUAL_BUILD_FLAGS=\"cmake_${CMAKE_CXX_COMPILER_ID}_${ALPHADUAL_FLAGS_DESC}\"")

# Put the executable next to the scripts, whatever the generator/configuration.
set_target_properties(alphadual PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
foreach(cfg Debug Release RelWithDebInfo MinSizeRel)
  string(TOUPPER ${cfg} CFG)
  set_target_properties(alphadual PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY_${CFG} ${CMAKE_CURRENT_SOURCE_DIR})
endforeach()
