cmake_minimum_required(VERSION 3.6)
project(FPMv1)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories(include)
set(SOURCE_FILES src/main.cpp
        src/pre.cpp
            src/read_config.cpp
            src/read_sgrid.cpp
            src/create_panel.cpp
            src/create_wake.cpp
            src/create_field.cpp
            src/map_field.cpp
        src/solver.cpp
            #src/id_subpanel.cpp
            src/build_AIC.cpp
                src/infcBB.cpp
                #src/split_panel.cpp
                    src/infcBF.cpp
        src/infcF.cpp
            src/solve_body.cpp
            src/solve_field.cpp
        src/compute_fVars.cpp
                    #interp_ctv
                    #interp_sp
                        #interp
        src/compute_sVars.cpp
        src/post.cpp
            src/write_sp.cpp
        src/write_fv.cpp)
add_executable(frame ${SOURCE_FILES})
add_executable(celia ${SOURCE_FILES})