cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(Paco)

set(PACO_SRCS
Graph.cpp
Community.cpp
Surprise.cpp
AsymptoticSurprise.cpp
igraph_utils.cpp
SurpriseFunction.cpp
ConditionalSurpriseFunction.cpp
AsymptoticSurpriseFunction.cpp
SignificanceFunction.cpp
ModularityFunction.cpp
AsymptoticModularityFunction.cpp
RandomOptimizer.cpp
AnnealOptimizer.cpp
PartitionHelper.cpp
AgglomerativeOptimizer.cpp
KLDivergence.cpp
Timer.cpp
)

set(PACO_HDRS
Graph.h
Community.h
Surprise.h
AsymptoticSurprise.h
igraph_utils.h
QualityFunction.h
SurpriseFunction.h
ConditionalSurpriseFunction.h
AsymptoticSurprise.h
AsymptoticSurpriseFunction.h
SignificanceFunction.h
ModularityFunction.h
AsymptoticModularityFunction.h
QualityOptimizer.h
RandomOptimizer.h
AnnealOptimizer.h
PartitionHelper.h
AgglomerativeOptimizer.h
KLDivergence.h
Timer.h
)

add_library(PACO SHARED ${PACO_SRCS} ${PACO_HDRS})
target_link_libraries(PACO ${IGRAPH_LIBRARIES})

if(MATLAB_SUPPORT)
    add_mex(paco_mx paco.cpp)
    target_link_libraries(paco_mx PACO ${MATLAB_LIBRARIES})
endif(MATLAB_SUPPORT)

if(OCTAVE_SUPPORT)
octave_add_oct(paco_oct SOURCES paco.cpp LINK_LIBRARIES PACO EXTENSION "mex")
endif(OCTAVE_SUPPORT)

if(PYTHON_SUPPORT)
    set_source_files_properties(pypaco.pyx PROPERTIES CYTHON_IS_CXX 1 )
    cython_add_module(pypaco pypaco.pyx)
    target_link_libraries(pypaco PACO)
endif(PYTHON_SUPPORT)

add_executable(test_graphc test_graphc.cpp)
target_link_libraries(test_graphc PACO)

add_executable(test_paco_multiple test_paco_multiple.cpp)
target_link_libraries(test_paco_multiple PACO)

add_executable(test_igraph_matrix_view test_igraph_matrix_view.cpp)
target_link_libraries(test_igraph_matrix_view PACO)

add_executable(test_igraph_matrix test_igraph_matrix.cpp)
target_link_libraries(test_igraph_matrix PACO)

add_executable(test_quality_function test_quality_function.cpp)
target_link_libraries(test_quality_function PACO)

add_executable(test_partition_helper test_partition_helper.cpp)
target_link_libraries(test_partition_helper PACO)

add_executable(test_optimizer test_optimizer.cpp)
target_link_libraries(test_optimizer PACO)

add_executable(test_hypergeometric test_hypergeometric.cpp)
target_link_libraries(test_hypergeometric PACO)