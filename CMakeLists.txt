# (C) British Crown Copyright 2024 Met Office
############################################################################################
# nemo-feedback

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

################################################################################
# Initialise project nemo-feedback

project( nemo-feedback VERSION 0.4.0 LANGUAGES CXX Fortran )
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package( ecbuild QUIET )
include( ecbuild_system NO_POLICY_SCOPE )
ecbuild_declare_project()

set( CMAKE_DIRECTORY_LABELS "nemo-feedback" )

################################################################################
# Required packages

find_package( eckit 1.16 REQUIRED )
ecbuild_debug( "   eckit_FEATURES : [${eckit_FEATURES}]" )

find_package( oops REQUIRED )
ecbuild_debug( "   oops_FEATURES : [${oops_FEATURES}]" )

find_package( ioda REQUIRED )
ecbuild_debug( "   ioda_FEATURES : [${ioda_FEATURES}]" )

find_package( ufo REQUIRED )
ecbuild_debug( "   ufo_FEATURES : [${ufo_FEATURES}]" )

ecbuild_find_package( NAME NetCDF COMPONENTS CXX Fortran )
find_package( NetCDF COMPONENTS CXX)
ecbuild_debug( "   NetCDF_FEATURES: [${NetCDF_FEATURES}]" )


################################################################################
# Features that can be enabled / disabled with -DENABLE_<FEATURE>

#include( features/BOUNDSCHECKING )

################################################################################
# sources

include_directories ( ${CMAKE_CURRENT_SOURCE_DIR}/src
                      ${CMAKE_CURRENT_BINARY_DIR}/src
                    )

add_subdirectory( src )

################################################################################
# Export and summarize

ecbuild_install_project( NAME nemo_feedback )

ecbuild_print_summary()

