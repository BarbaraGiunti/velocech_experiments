cmake_minimum_required(VERSION 3.16)

project(alpha_benchmark LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# GUDHI source checkout
set(GUDHI_DIR "$ENV{HOME}/projects/gudhi-devel")

# GUDHI and Eigen headers
include_directories(
    "${GUDHI_DIR}/src"
    "${GUDHI_DIR}/src/Alpha_complex/include"
    "${GUDHI_DIR}/src/Simplex_tree/include"
    "C:/msys64/ucrt64/include/eigen3"
    "C:/msys64/home/bg791928/projects/gudhi-devel/src/common/include"
)

add_executable(alpha_benchmark
    alpha_benchmark.cpp
)


target_include_directories(alpha_benchmark PRIVATE
    C:/msys64/home/bg791928/projects/gudhi-devel/src
    C:/msys64/home/bg791928/projects/gudhi-devel/src/Alpha_complex/include
    C:/msys64/home/bg791928/projects/gudhi-devel/src/Simplex_tree/include
    C:/msys64/home/bg791928/projects/gudhi-devel/src/common/include
    C:/msys64/ucrt64/include/eigen3
)

target_link_libraries(alpha_benchmark PRIVATE
    gmp
    mpfr
)
