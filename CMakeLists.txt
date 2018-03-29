cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(Paco)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
set(PROJECT_BINARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/build/")
set(LIBRARY_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/build/" CACHE PATH "Single directory for all libraries")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/build/" CACHE PATH "Single directory for all executables")

# Check for the support of C++11 compliant compilers
include(CheckCXXCompilerFlag)

if (UNIX)
set(CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
endif(UNIX)
# Specify if matlab support and python wrappers are to be compiled. Default OFF
# MATLAB support needs MEX (Tested on Ubuntu 14.04, Matlab R2013a)
option(MATLAB_SUPPORT "Mex with MATLAB" OFF)
option(OCTAVE_SUPPORT "Oct with Octave mkoctfile" OFF)
option(SAMPLE_LANDSCAPE "Enable the SAMPLE_LANDSCAPE option, to write down at every stage of the optimization process, the quality function and the membership")
option(EXPERIMENTAL_FEATURES "Enable compilation of some experimental features, default FALSE" OFF)
option(COMPILE_TESTS "Compile all debugging tests" OFF)

if(MATLAB_SUPPORT)
	if (WIN32)
		#set(MATLAB_ADDITIONAL_VERSIONS "R2016a=9.0.1")
		#set(MATLAB_ROOT "C:\\Program\ Files\\MATLAB\\R2016a")
		message(STATUS ${MATLAB_ROOT})
		find_package(MyMatlab REQUIRED)
		message(STATUS "Matlab include dir= ${MATLAB_INCLUDE_DIR}")
		message(STATUS "Matlab MX library = ${_MATLABLIB}")
		add_definitions(/DMATLAB_MEX_FILE) #define matlab macros
		add_definitions(/DMX_COMPAT_64)
		message(STATUS ${MATLAB_INCLUDE_DIRS})
		message(STATUS ${MATLAB_LIBRARIES})
		include_directories(${MATLAB_INCLUDE_DIRS})
	elseif(UNIX)
		find_package(MyMatlab)
		message(STATUS "Matlab include directories: ${MATLAB_INCLUDE_DIRS}")
	endif(WIN32)
	
    
    include_directories(${MATLAB_INCLUDE_DIRS})
    if (MATLAB_NOT_FOUND)
        message(FATAL "Error matlab not found")
    endif()
    add_definitions("-DMATLAB_SUPPORT")
endif()

if(OCTAVE_SUPPORT)
    find_package(OCTAVE)
    include_directories(${OCTAVE_INCLUDE_DIRS})
    add_definitions("-DOCTAVE_SUPPORT")
endif(OCTAVE_SUPPORT)


# PYTHON_SUPPORT needs Cython (install it via "pip install cython")
option(PYTHON_SUPPORT "Python support through Cython extensions" OFF)
if(PYTHON_SUPPORT)
    set(Python_ADDITIONAL_VERSIONS "2.7.12")
    find_package(PythonInterp)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    find_package(PythonLibs 2.7.12 EXACT REQUIRED)
    find_package(Numpy)
    include_directories(${PYTHON_NUMPY_INCLUDE_DIR})
        # Print some informations on python libraries linked and included
        message(STATUS "Python libraries are ${PYTHON_LIBRARIES}")
        message(STATUS "Python libraries are included in ${PYTHON_INCLUDE_DIRS}")
        message(STATUS "Python versions found is ${PYTHONLIBS_VERSION_STRING}")
        include_directories("${PYTHON_INCLUDE_DIR}")
    include(UseCython)
    include_directories("${PYTHON_INCLUDE_DIR}/numpy")
    message(STATUS "${PYTHON_INCLUDE_DIR}/numpy")
#    add_custom_target( ReplicatePythonSourceTree ALL ${CMAKE_COMMAND} -P
#      ${CMAKE_CURRENT_SOURCE_DIR}/cmake/ReplicatePythonSourceTree.cmake
#      ${CMAKE_CURRENT_BINARY_DIR}
#      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    # This definition avoids a stupid warning message
    #add_definitions("-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION")
endif(PYTHON_SUPPORT)


if(SAMPLE_LANDSCAPE)
    add_definitions("-DSAMPLE_LANDSCAPE")
endif()

if(EXPERIMENTAL_FEATURES)
    add_definitions("-DEXPERIMENTAL_FEATURES")
endif()

if(PYTHON_SUPPORT)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
    find_package(PythonLibs REQUIRED )
        # Print some informations on python libraries linked and included
        message(STATUS "Python libraries are ${PYTHON_LIBRARIES}")
        message(STATUS "Python libraries are included in ${PYTHON_INCLUDE_DIRS}")
        message(STATUS "Python versions found is ${PYTHONLIBS_VERSION_STRING}")
        include_directories("${PYTHON_INCLUDE_DIR}")
    include(UseCython)
    include_directories()
    include_directories("${PYTHON_INCLUDE_DIR}/numpy")
    message(STATUS "${PYTHON_INCLUDE_DIR}/numpy")
    add_custom_target( ReplicatePythonSourceTree ALL ${CMAKE_COMMAND} -P
      ${CMAKE_CURRENT_SOURCE_DIR}/cmake/ReplicatePythonSourceTree.cmake
      ${CMAKE_CURRENT_BINARY_DIR}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    # This definition avoids a stupid warning message
    #add_definitions("-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION")
endif(PYTHON_SUPPORT)

if (MSVC)
set(IGRAPH_INCLUDES "${PROJECT_SOURCE_DIR}/igraph-0.7.1-msvc/igraph-0.7.1-msvc/include")
set(IGRAPH_LIBRARIES "${PROJECT_SOURCE_DIR}/igraph-0.7.1-msvc/igraph-0.7.1-msvc/Release/igraph.lib")
else (MSVC)
find_package(IGraph)
endif(MSVC)

include_directories(${IGRAPH_INCLUDES})

# Add the subdirectory with the source code
add_subdirectory(src)
