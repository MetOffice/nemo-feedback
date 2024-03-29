# (C) British Crown Copyright 2024 Met Office

cmake_minimum_required( VERSION 3.18 FATAL_ERROR )

project( nemo-feedback-ci VERSION 0.1 LANGUAGES C CXX Fortran )

set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )
set( ENABLE_OMP ON CACHE BOOL "Compile with OpenMP" )

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/jedicmake" EXCLUDE_FROM_ALL)
if(NOT DEFINED jedicmake_DIR)
  set(jedicmake_DIR "${CMAKE_BINARY_DIR}/jedicmake")
endif()

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/oops" EXCLUDE_FROM_ALL)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/ioda" EXCLUDE_FROM_ALL)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/ufo" EXCLUDE_FROM_ALL)
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/nemo-feedback")
enable_testing()
