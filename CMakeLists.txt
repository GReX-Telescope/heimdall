cmake_minimum_required(VERSION 3.18)

# Setup project definition
project(heimdall VERSION 0.1 LANGUAGES CUDA CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# CUDA Configuration
set(CUDA_NVCC_FLAGS -O3 -use_fast_math -lineinfo -src-in-ptx)

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 75)
endif()

# Define dependencies
include(GNUInstallDirs)
# CUDAToolkit
find_package(CUDAToolkit REQUIRED)
# Thrust
find_package(Thrust REQUIRED CONFIG)
thrust_create_target(Thrust)
# Boost
find_package(Boost REQUIRED)
# Dedisp
find_package(Dedisp REQUIRED)
# PSRDada

# Add CUDA headers for non cu/cuh
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}) 

# Add subdirectories
add_subdirectory(Applications)
add_subdirectory(Formats)
add_subdirectory(Network)
add_subdirectory(Pipeline)