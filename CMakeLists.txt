cmake_minimum_required(VERSION 3.6)
project(FPMv1)

#IF(UNIX AND NOT APPLE)
#IF(WIN32)
IF(UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DON_UNIX")
ENDIF()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

IF (UNIX)
    include_directories(include)
ELSEIF(WIN32)
    include_directories(include; C:/eigen-3.3.3/Eigen/)
ENDIF()
set(SOURCE_FILES src/main.cpp
        src/pre.cpp
        src/read_config.cpp
        src/read_sgrid.cpp
        src/create_panel.cpp
        src/create_wake.cpp
        src/create_field.cpp
        src/map_field.cpp
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

IF(CMAKE_CXX_FLAGS_RELEASE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCMAKE_CXX_FLAGS=-O3")
    add_executable(celia ${SOURCE_FILES})
ELSEIF(CMAKE_CXX_FLAGS_DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCMAKE_CXX_FLAGS=-G -DCMAKE_CXX_FLAGS=-O2")
    add_executable(frame ${SOURCE_FILES})
ENDIF()

