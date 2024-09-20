cmake_minimum_required(VERSION 3.29)
project(PoissonDiskSamplingCPP)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add optimization flag for release builds
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

add_executable(PoissonDiskSamplingCPP main.cpp
        PoissonDiskSampling.cpp
        PoissonDiskSampling.h)
