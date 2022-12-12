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

target_link_libraries(hdpipeline hdnetwork hdformats Thrust ${Boost_LIBRARIES} ${DEDISP_LIBRARIES})

target_include_directories(hdpipeline PRIVATE ${DEDISP_INCLUDE_DIR})
target_include_directories(hdpipeline PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(hdpipeline PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../Network")
target_include_directories(hdpipeline PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/../Formats")

set_target_properties(hdpipeline PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

target_compile_options(hdpipeline PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--extended-lambda>)