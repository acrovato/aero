# FindEIGEN.cmake - try to find Eigen headers
# ----------------------------------------------------------------------------   
# output:
#    EIGEN_FOUND             : TRUE/FALSE
#    EIGEN_INCLUDE_DIRS      : where the Eigen/*.h are              [cached]
# ----------------------------------------------------------------------------
# autodetection:
#    set "CMAKE_INCLUDE_PATH=path/to/eigen/"
#     or "INCLUDE=path/to/eigen/"
# ----------------------------------------------------------------------------   

find_path(EIGEN_INCLUDE_DIRS "Eigen/Dense" PATHS "/usr/include/eigen3")

# handle the QUIETLY and REQUIRED arguments and set EIGEN_FOUND to TRUE
# if all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(EIGEN DEFAULT_MSG 
				                  EIGEN_INCLUDE_DIRS)

#MESSAGE(STATUS "EIGEN_FOUND = ${EIGEN_FOUND}")
