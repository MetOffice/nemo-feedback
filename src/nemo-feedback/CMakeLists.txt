# (C) Crown Copyright 2021, the Met Office. All rights reserved.
#
# Refer to COPYRIGHT.txt of this distribution for details.

list( APPEND nemo_fdbk_src_files
    instantiateObsFilterFactory.h
    NemoFeedback.cc
    NemoFeedback.h
    NemoFeedbackParameterTraitsOutputDtype.h
    NemoFeedbackParameterTraitsOutputDtype.cc
    NemoFeedbackParameters.h
    NemoFeedbackDataCreator.h
    NemoFeedbackDataCreator.cc
    feedback_io/DataIndexer.h
    feedback_io/DataIndexer.cc
    feedback_io/Data.h
    feedback_io/Data.cc
    feedback_io/Utils.h
    feedback_io/Utils.cc
    feedback_io/Writer.h
    feedback_io/Writer.cc
)

ecbuild_add_library( TARGET nemo_feedback
                     SOURCES ${nemo_fdbk_src_files}
                     INSTALL_HEADERS LISTED
                     HEADER_DESTINATION ${INSTALL_INCLUDE_DIR}/nemo-feedback
                   )
set_target_properties(nemo_feedback PROPERTIES CXX_STANDARD 17)
target_link_libraries( nemo_feedback PUBLIC ufo ioda oops NetCDF::NetCDF_CXX)
target_include_directories(nemo_feedback PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>
                                                $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
