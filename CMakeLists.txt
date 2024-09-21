cmake_minimum_required(VERSION 3.29)
project(PoissonDiskSamplingCPP)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add optimization flag for release builds
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

find_package(Boost REQUIRED COMPONENTS system filesystem)

add_executable(PoissonDiskSamplingCPP main.cpp
        poisson_disk_sampling.cpp
        poisson_disk_sampling.h
        distance_constrained_sampling.cpp
        distance_constrained_sampling.h
        distance_metrics.h
        sampling_interface.h)

target_link_libraries(PoissonDiskSamplingCPP PRIVATE Boost::system Boost::filesystem)