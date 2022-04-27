# (C) Copyright 2022 Met Office
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.


cmake_minimum_required( VERSION 3.18 FATAL_ERROR )

project( nemo-feedback-ci VERSION 0.1 LANGUAGES C CXX Fortran )

find_package(ecbuild)
set( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake;${CMAKE_CURRENT_SOURCE_DIR}/ecbuild/cmake;${CMAKE_MODULE_PATH}")

include( ecbuild_bundle )

set( ENABLE_MPI ON CACHE BOOL "Compile with MPI" )
set( ENABLE_OMP ON CACHE BOOL "Compile with OpenMP" )
# Prevents CMake from being re-run after each commit (useful for development,
# not necessarily for production).
set( ECBUILD_RECORD_GIT_COMMIT_SHA1 OFF CACHE BOOL
     "When building ecbuild projects that are Git repos, create variables recording the full and short Git revision" )

ecbuild_bundle_initialize()

ecbuild_requires_macro_version( 2.7 )

ecbuild_bundle( PROJECT jedicmake SOURCE $ENV{CMAKE_CURRENT_SOURCE_DIR}/jedicmake )
include( jedicmake/cmake/Functions/git_functions.cmake )

ecbuild_bundle( PROJECT eckit SOURCE $ENV{CMAKE_CURRENT_SOURCE_DIR}/eckit )
ecbuild_bundle( PROJECT fckit SOURCE $ENV{CMAKE_CURRENT_SOURCE_DIR}/fckit )
ecbuild_bundle( PROJECT atlas SOURCE $ENV{CMAKE_CURRENT_SOURCE_DIR}/atlas )
ecbuild_bundle( PROJECT oops  SOURCE $ENV{CMAKE_CURRENT_SOURCE_DIR}/oops )
ecbuild_bundle( PROJECT ioda  SOURCE $ENV{CMAKE_CURRENT_SOURCE_DIR}/ioda )
ecbuild_bundle( PROJECT ufo  SOURCE $ENV{CMAKE_CURRENT_SOURCE_DIR}/ufo )
ecbuild_bundle( PROJECT atlas-orca SOURCE $ENV{CMAKE_CURRENT_SOURCE_DIR}/atlas-orca )
ecbuild_bundle( PROJECT nemo-feedback SOURCE /usr/local/src/nemo-feedback )
ecbuild_bundle( PROJECT orca-jedi SOURCE $ENV{CMAKE_CURRENT_SOURCE_DIR}/orca-jedi )

ecbuild_bundle_finalize()
