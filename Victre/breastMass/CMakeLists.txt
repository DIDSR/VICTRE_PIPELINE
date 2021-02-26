cmake_minimum_required(VERSION 2.8)

PROJECT(breastMassProject)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(OpenMP)
if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories("$(PROJECT_SOURCE_DIR)")

SET(CMAKE_BUILD_TYPE "Release")
SET(CMAKE_CXX_FLAGS  "-std=c++0x -fext-numeric-literals ${CMAKE_CXX_FLAGS}")

add_executable(breastMass breastMass.cxx)
LINK_DIRECTORIES("$(LAPACK_LIB)")
LINK_DIRECTORIES("$(BOOST_OPTIONS_LIB)")
target_link_libraries(breastMass lapack blas boost_program_options ${VTK_LIBRARIES})
