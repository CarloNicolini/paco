cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(Paco)

# Include Eigen matrix libraries
include_directories(eigen)

add_subdirectory(paco)
add_subdirectory(lfrweighted)
