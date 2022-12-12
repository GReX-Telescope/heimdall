cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

# Setup project definition
project(heimdall VERSION 0.1 LANGUAGES CUDA CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CUDA_ARCHITECTURES 75 86)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Define dependencies
include(GNUInstallDirs)
# CUDAToolkit
find_package(CUDAToolkit REQUIRED)
# Thrust
find_package(Thrust REQUIRED CONFIG)
thrust_create_target(Thrust)
# Boost
find_package(Boost REQUIRED COMPONENTS system)
# Dedisp
add_subdirectory(dedisp)
set(DEDISP_LIBRARIES "${DEDISP_tdd_LIBRARY}")
# PSRDada
find_package(PSRDada)

if (PSRDada_FOUND)
  set(HAVE_PSRDADA on)
  add_compile_definitions(HAVE_PSRDADA)
  find_package(rdmacm REQUIRED)
  set(HAVE_RDMACM on)
  message(STATUS "Found PSRDada")
endif()

# Add headers
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}) 

# Add subdirectories
add_subdirectory(Applications)
add_subdirectory(Formats)
add_subdirectory(Network)
add_subdirectory(Pipeline)