cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(MPASMap)

set(CMAKE_CXX_STANDARD 11)

# Add .lib files
link_directories(
    "${CMAKE_SOURCE_DIR}/lib"
)
	
# Add header files
file(GLOB_RECURSE HEADER_FILES 
	${CMAKE_SOURCE_DIR}/src/*.h
	${CMAKE_SOURCE_DIR}/src/*.hpp)
	
# We need a CMAKE_DIR with some code to find external dependencies
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

#######################################
# LOOK for the packages that we need! #
#######################################

# OpenGL EGL
find_package(OpenGL REQUIRED COMPONENTS EGL)

# GLM
find_package(GLM REQUIRED)
message(STATUS "GLM included at ${GLM_INCLUDE_DIR}")

# GLAD
add_library(glad OBJECT third_party/glad/src/glad.c)
target_include_directories(glad PUBLIC third_party/glad/include/)
target_link_libraries(glad PUBLIC ${CMAKE_DL_LIBS})

add_executable(eglinfo src/eglinfo.c)
target_link_libraries(
    eglinfo
    PRIVATE
    OpenGL::EGL
    glad
)

# NETCDF
link_directories(/usr/local/netcdf/intel/17.0/mvapich2/2.3/4.6.1/lib)
link_libraries(netcdf)

# HDF5
link_directories(/usr/local/hdf5/intel/17.0/mvapich2/2.3/1.10.2/lib)
link_libraries(hdf5)

# MPI 
link_directories(/opt/mvapich2/intel/17.0/2.3.2/lib)
link_libraries(mpi)

# Define the include DIRs
include_directories(
	"${CMAKE_SOURCE_DIR}/src"
	"${CMAKE_SOURCE_DIR}/include"
	"/usr/local/netcdf/intel/17.0/mvapich2/2.3/4.6.1/include"
	"/usr/local/hdf5/intel/17.0/mvapich2/2.3/1.10.2/include"
	"/opt/mvapich2/intel/17.0/2.3.2/include"
)

# Define the executable
add_executable(${PROJECT_NAME} ${HEADER_FILES} src/main.cpp)

# Put all libraries into a variable
set(LIBS OpenGL glad GL X11 pthread Xrandr Xi dl netcdf hdf5 mpi OpenGL::EGL)

# Define the link libraries
target_link_libraries(${PROJECT_NAME} ${LIBS})
