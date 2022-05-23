cmake_minimum_required(VERSION 3.18)

# Setup project definition
project(heimdall VERSION 0.1 LANGUAGES CUDA CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES "35;50;75")
endif()

if(NOT DEFINED CMAKE_CUDA_STANDARD)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
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
find_package(PSRDada)

if (PSRDADA_FOUND)
  set(HAVE_PSRDADA on)
  add_compile_definitions(HAVE_PSRDADA)
endif()


# Add headers
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}) 

# Add subdirectories
add_subdirectory(Applications)
add_subdirectory(Formats)
add_subdirectory(Network)
add_subdirectory(Pipeline)