
############################################################################################
# nemo-feedback

cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

find_package( ecbuild 3.4 REQUIRED HINTS ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild )
include( ecbuild_system NO_POLICY_SCOPE )

################################################################################
# Initialise project nemo-feedback

project( nemo-feedback LANGUAGES CXX )
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

find_package( orca-jedi )

#ecbuild_find_package( NAME NetCDF COMPONENTS CXX)
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

