set(format_sources SigprocFile.cpp)

if (PSRDada_FOUND)
    set(format_sources ${format_sources} PSRDadaRingBuffer.cpp)
endif()

add_library(hdformats ${format_sources})
target_include_directories(hdformats PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

if (PSRDada_FOUND)
    target_link_libraries(hdformats psrdada)
endif()