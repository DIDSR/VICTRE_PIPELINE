cmake_minimum_required(VERSION 2.8)

PROJECT(breastPhantomProject)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(OpenMP)
if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

#add_subdirectory(vtkPoissonReconstruction)

#include_directories("/opt/local/include")
#include_directories("/usr/local/include")
include_directories("$(PROJECT_SOURCE_DIR)")
#include_directories("$(PROJECT_SOURCE_DIR)/vtkPoissonReconstruction")

add_library(perlinNoise perlinNoise.cxx)
add_library(duct duct.cxx)
add_library(artery artery.cxx)
add_library(vein vein.cxx)
add_library(createDuct createDuct.cxx)
add_library(createArtery createArtery.cxx)
add_library(createVein createVein.cxx)

SET(CMAKE_BUILD_TYPE "Release")
#SET(DISABLE_IRRITATING_WARNING "-Wno-unnamed-type-template-args")
SET(CMAKE_CXX_FLAGS  "-std=c++0x -fext-numeric-literals ${CMAKE_CXX_FLAGS}")

#LINK_DIRECTORIES("/opt/local/lib")
#LINK_DIRECTORIES("$(PROJECT_SOURCE_DIR)/vtkPoissonReconstruction")
LINK_DIRECTORIES("$(LAPACK_LIB)")
LINK_DIRECTORIES("$(BOOST_OPTIONS_LIB)")
#add_executable(breastPhantomMain MACOSX_BUNDLE breastPhantom.cxx)
add_executable(breastPhantomMain breastPhantom.cxx)

target_link_libraries(breastPhantomMain perlinNoise createDuct createArtery createVein duct artery vein z lapack blas boost_program_options ${VTK_LIBRARIES})
#target_link_libraries(breastPhantomMain PoissonReconstruction createDuct createArtery createVein duct artery vein z lapack blas boost_program_options-mt ${VTK_LIBRARIES})
