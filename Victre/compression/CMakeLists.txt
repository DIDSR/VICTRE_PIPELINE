cmake_minimum_required(VERSION 2.8)

PROJECT(breastCompressProject)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(OpenMP)
if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories("$(PROJECT_SOURCE_DIR)")

SET(CMAKE_BUILD_TYPE "Release")
SET(CMAKE_CXX_FLAGS  "-std=c++0x ${CMAKE_CXX_FLAGS}")

add_executable(breastCompress breastCompress.cxx)
LINK_DIRECTORIES("~/libraries/lapack-master/build/lib")
LINK_DIRECTORIES("/home/kenny.cha/boost_1_66_0/bin.v2/libs/program_options/build/gcc-gnu-4.8.5/release/link-static/threadapi-pthread/threading-multi/")
target_link_libraries(breastCompress z tet lapack blas boost_program_options ${VTK_LIBRARIES})
