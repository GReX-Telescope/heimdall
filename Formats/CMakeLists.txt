add_library(hdformats SigprocFile.cpp PSRDadaRingBuffer.cpp)
target_include_directories(hdformats PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})