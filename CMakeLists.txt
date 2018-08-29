cmake_minimum_required(VERSION 3.6)
project(FPMv1)

# Check OS
#IF(UNIX AND NOT APPLE)
#IF(WIN32)
IF(UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DON_UNIX")
ENDIF()

# Check for C++11
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
    message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

# Add a automated way to find Eigen
# Set path for header and source files
IF (UNIX)
    include_directories(include; /usr/local/include/eigen-3.3.5/)
ELSEIF(WIN32)
    include_directories(include; C:/eigen-3.3.5/)
ENDIF()
set(SOURCE_FILES src/main.cpp
        src/pre.cpp
        src/read_config.cpp
        src/read_sgrid.cpp
        src/create_panel.cpp
        src/create_wake.cpp
        src/create_field.cpp
        src/map_field.cpp
        src/map_derivatives.cpp
        src/cast_ray_pip.cpp
        src/solver.cpp
        src/id_subpanel.cpp
        src/build_AIC.cpp
        src/infcB.cpp
        src/split_panel.cpp
        src/infcFF.cpp
        src/infcFB.cpp
        src/solve_body.cpp
        src/solve_field.cpp
        src/compute_fVars.cpp
        src/interp_ctv
        src/interp_sp
        src/interp
        src/compute_sVars.cpp
        src/post.cpp
        src/write_sp.cpp
        src/write_fv.cpp)

# Specifiy default build type
if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release CACHE STRING
       "Choose the type of build, options are: None Debug Release RelWithDebInfo
MinSizeRel."
       FORCE )
endif()
message(STATUS "Build type: " ${CMAKE_BUILD_TYPE})

# Set verbose mode
if(CMAKE_BUILD_TYPE MATCHES Debug OR CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DVERBOSE")
endif()

# Specifiy exe location
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# Set Release or Debug option
IF(CMAKE_CXX_FLAGS_RELEASE OR CMAKE_CXX_FLAGS_RELWITHDEBINFO)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCMAKE_CXX_FLAGS=-O3")
    add_executable(celia ${SOURCE_FILES})
ELSEIF(CMAKE_CXX_FLAGS_DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCMAKE_CXX_FLAGS=-G -DCMAKE_CXX_FLAGS=-O2")
    add_executable(frame ${SOURCE_FILES})
ENDIF()

