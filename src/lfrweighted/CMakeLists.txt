cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(LFRW)

set(LFRW_SRCS
benchm.cpp
binary_benchm.cpp
cast.cpp
cc.cpp
combinatorics.cpp
histograms.cpp
print.cpp
random.cpp
set_parameters.cpp
)

set(LFRW_HDRS
benchm.h
binary_benchm.h
cast.h
cc.h
combinatorics.h
histograms.h
print.h
random.h
set_parameters.h
standard_include.h
)


add_library(LFRW SHARED ${LFRW_SRCS} ${LFRW_HDRS})
target_link_libraries(LFRW ${IGRAPH_LIBRARIES})

add_executable(test_benchm test_benchm.cpp)
target_link_libraries(test_benchm LFRW)

# if(MATLAB_SUPPORT)
#     # To avoid a stupid linking problem with incompatible libstdc++ versions
#     # http://stackoverflow.com/questions/19386651/how-to-fix-usr-lib-libstdc-so-6-version-glibcxx-3-4-15-not-found
#     #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -static-libstdc++")
#     add_mex(lfrw_mx lfrw.cpp)
#     target_link_libraries(lfrw_mx LFRW ${MATLAB_LIBRARIES})
# endif(MATLAB_SUPPORT)
