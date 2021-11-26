# (C) Crown Copyright 2021, the Met Office. All rights reserved.
#
# Refer to COPYRIGHT.txt of this distribution for details.

list( APPEND nemo_fdbk_src_files
    instantiateObsFilterFactory.h
    NemoFeedback.cc
    NemoFeedback.h
    NemoFeedbackParameters.h
    NemoFeedbackWriter.h
    NemoFeedbackWriter.cc
)

ecbuild_add_library( TARGET nemo_feedback
                     SOURCES ${nemo_fdbk_src_files}
                     INSTALL_HEADERS LISTED
                     HEADER_DESTINATION ${INSTALL_INCLUDE_DIR}/nemo-feedback
                   )
target_link_libraries( nemo_feedback PUBLIC ufo ioda oops NetCDF::NetCDF_CXX)
target_include_directories(nemo_feedback PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
                                                $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
