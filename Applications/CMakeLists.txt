
# List all the executables we'll create
add_executable(heimdall heimdall.cpp)
add_executable(coincidencer coincidencer.cpp Candidates.cpp)
add_executable(coincidencer_client coincidencer_client.cpp)
add_executable(candidate_profiler candidate_profiler.cpp)
add_executable(fil2pgm fil2pgm.cpp)
add_executable(generate_dmlist generate_dmlist.cpp)

# heimdall - main executable
target_include_directories(heimdall PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../Pipeline")
target_include_directories(heimdall PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../Formats")
target_link_libraries(heimdall hdpipeline hdnetwork hdformats)
if (PSRDADA_FOUND)
    target_link_libraries(heimdall psrdada)
endif()
install(TARGETS heimdall DESTINATION bin)

# candidate_profiler
target_link_libraries(candidate_profiler hdformats dedisp)
target_include_directories(candidate_profiler PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../Formats")
install(TARGETS candidate_profiler DESTINATION bin)

# coincidencer_client
target_link_libraries(coincidencer_client hdnetwork)
target_include_directories(coincidencer_client PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../Network")
install(TARGETS coincidencer_client DESTINATION bin)

# coincidencer
target_link_libraries(coincidencer hdnetwork)
target_include_directories(coincidencer PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../Network")
install(TARGETS coincidencer DESTINATION bin)

# fil2pgm
target_link_libraries(fil2pgm hdformats)
target_include_directories(fil2pgm PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../Formats")
install(TARGETS fil2pgm DESTINATION bin)

# generate_dmlist
target_link_libraries(generate_dmlist hdformats hdpipeline dedisp)
target_include_directories(generate_dmlist PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../Formats")
target_include_directories(generate_dmlist PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../Pipeline")
install(TARGETS generate_dmlist DESTINATION bin)