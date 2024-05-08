cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

# Setup project definition
project(heimdall VERSION 0.1 LANGUAGES CXX CUDA)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "-DBOOST_LOG_DYN_LINK")

set(CMAKE_CUDA_ARCHITECTURES 60 61 62 70 72 75 86)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Define dependencies
include(GNUInstallDirs)
# CUDA
find_package(CUDAToolkit REQUIRED)
# Boost
find_package(Boost REQUIRED COMPONENTS log log_setup REQUIRED)
# Dedisp
add_subdirectory(dedisp)
# PSRDada (with pkg-config)
find_package(PkgConfig REQUIRED)
pkg_search_module(PSRDada psrdada)
if (PSRDada_FOUND)
  set(HAVE_PSRDADA on)
  add_compile_definitions(HAVE_PSRDADA)
  message(STATUS "Found PSRDada")
endif()

# Add headers
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}) 

# Add subdirectories
add_subdirectory(Applications)
add_subdirectory(Formats)
add_subdirectory(Network)
add_subdirectory(Pipeline)
