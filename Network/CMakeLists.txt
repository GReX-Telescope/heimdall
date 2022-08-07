add_library(hdnetwork ClientSocket.cpp Socket.cpp ServerSocket.cpp)
target_include_directories(hdnetwork PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})