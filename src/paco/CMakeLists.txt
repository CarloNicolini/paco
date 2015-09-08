cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(Paco)

set(PACO_SRCS
Graph.cpp
Community.cpp
Surprise.cpp
igraph_utils.cpp
SurpriseFunction.cpp
AsymptoticSurprise.cpp
AsymptoticSurpriseFunction.cpp
ModularityFunction.cpp
SignificanceFunction.cpp
RandomOptimizer.cpp
AnnealOptimizer.cpp
PartitionHelper.cpp
AgglomerativeOptimizer.cpp
KLDivergence.cpp
)

set(PACO_HDRS
Graph.h
Community.h
Surprise.h
igraph_utils.h
QualityFunction.h
SurpriseFunction.h
AsymptoticSurprise.h
AsymptoticSurpriseFunction.h
SignificanceFunction.h
ModularityFunction.h
QualityOptimizer.h
RandomOptimizer.h
AnnealOptimizer.h
PartitionHelper.h
AgglomerativeOptimizer.h
KLDivergence.h
)

add_library(PACO SHARED ${PACO_SRCS} ${PACO_HDRS})
target_link_libraries(PACO ${IGRAPH_LIBRARIES})

add_executable(test_graphc test_graphc.cpp)
target_link_libraries(test_graphc PACO)

add_executable(test_igraph_matrix test_igraph_matrix.cpp)
target_link_libraries(test_igraph_matrix PACO)

add_executable(test_quality_function test_quality_function.cpp)
target_link_libraries(test_quality_function PACO)

add_executable(test_partition_helper test_partition_helper.cpp)
target_link_libraries(test_partition_helper PACO)

add_executable(test_optimizer test_optimizer.cpp)
target_link_libraries(test_optimizer PACO)

if(MATLAB_SUPPORT)
    add_mex(paco paco.cpp)
    target_link_libraries(paco PACO ${MATLAB_LIBRARIES})
endif(MATLAB_SUPPORT)

if(OCTAVE_SUPPORT)
    octave_add_oct(paco_oct SOURCES paco.cpp LINK_LIBRARIES PACO EXTENSION "mex")
endif(OCTAVE_SUPPORT)
