
# List all the executables we'll create
add_executable(heimdall heimdall.cpp)
add_executable(coincidencer coincidencer.cpp Candidates.cpp)
add_executable(coincidencer_client coincidencer_client.cpp)
add_executable(candidate_profiler candidate_profiler.cpp)
add_executable(fil2pgm fil2pgm.cpp)
add_executable(generate_dmlist generate_dmlist.cpp)

# Heimdall main executable
target_include_directories(heimdall PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../Pipeline")
target_include_directories(heimdall PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../Formats")
target_link_libraries(heimdall STATIC hdformats)
target_link_libraries(heimdall STATIC hdnetwork)
target_link_libraries(heimdall STATIC hdpipeline)

target_include_directories(candidate_profiler PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../Formats")