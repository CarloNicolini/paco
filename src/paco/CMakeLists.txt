cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(Paco)

set(PACO_SRCS
Graph.cpp
Community.cpp
Surprise.cpp
AsymptoticSurprise.cpp
igraph_utils.cpp

SurpriseFunction.cpp
AsymptoticSurpriseFunction.cpp
SignificanceFunction.cpp
ModularityFunction.cpp

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
AsymptoticSurprise.h
AsymptoticSurpriseFunction.h
SignificanceFunction.h
ModularityFunction.cpp

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

add_executable(test_graphc test_graphc.cpp)
target_link_libraries(test_graphc PACO)

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

#if(MATLAB_SUPPORT)
# To avoid a stupid linking problem with incompatible libstdc++ versions
# http://stackoverflow.com/questions/19386651/how-to-fix-usr-lib-libstdc-so-6-version-glibcxx-3-4-15-not-found
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++")
#add_library(paco_mx SHARED paco.cpp ${PACO_SRCS})
#target_link_libraries(paco_mx ${MATLAB_LIBRARIES} ${IGRAPH_LIBRARIES})
#set_target_properties(paco_mx PROPERTIES SUFFIX ".${MATLAB_MEX_EXT}")
# Remove the "lib" prefix in front of mex files
#set_target_properties(paco_mx PROPERTIES PREFIX "")
#endif(MATLAB_SUPPORT)

if(OCTAVE_SUPPORT)
    octave_add_oct(paco_oct SOURCES paco.cpp LINK_LIBRARIES PACO EXTENSION "mex")
endif(OCTAVE_SUPPORT)

### EXPERIMENTAL FEATURES, TO ALLOW LINEAR COMBINATIONS OF QUALITY FUNCTIONS ##
if (EXPERIMENTAL)
        set(PACO_SRCS_EXPERIMENTAL
        ${PACO_SRCS}
        Graph.cpp
        Community.cpp
        AsymptoticSurprise.cpp
        KLDivergence.cpp
        Surprise.cpp
        igraph_utils.cpp
        SurpriseFunction2.cpp
        AsymptoticSurpriseFunction2.cpp
        )

        set(PACO_HDRS_EXPERIMENTAL
        ${PACO_HDRS}
        Graph.h
        Community.h
        Surprise.h
        igraph_utils.cpp
        QualityFunctionImpl.h
        SurpriseFunction2.h
        AsymptoticSurpriseFunction2.h
        UnionFunction.h
        SumFunction.h
        QualityFunction2.h
        )

        add_library(PACO_EXPERIMENTAL ${PACO_SRCS_EXPERIMENTAL} ${PACO_HDRS_EXPERIMENTAL})
        target_link_libraries(PACO_EXPERIMENTAL ${IGRAPH_LIBRARIES})
        add_executable(test_paco_exp test_paco_exp.cpp)
        target_link_libraries(test_paco_exp PACO_EXPERIMENTAL)
endif(EXPERIMENTAL)
###
