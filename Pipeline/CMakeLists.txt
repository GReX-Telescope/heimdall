add_library(hdpipeline
    default_params.cpp
    error.cpp
    parse_command_line.cpp
    clean_filterbank_rfi.cu
    get_rms.cu
    matched_filter.cu
    remove_baseline.cu
    find_giants.cu
    label_candidate_clusters.cu
    merge_candidates.cu
    pipeline.cu
    measure_bandpass.cu
    median_filter.cu
    matched_filter.cu)

target_link_libraries(hdpipeline Thrust Boost Dedisp)

target_include_directories(hdpipeline PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(hdpipeline PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../Network")
target_include_directories(hdpipeline PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../Formats")

target_compile_options(hdpipeline PRIVATE -Wall)