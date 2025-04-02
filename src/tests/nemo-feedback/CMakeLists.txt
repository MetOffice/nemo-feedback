# (C) British Crown Copyright 2024 Met Office

include_directories( ${PROJECT_SOURCE_DIR}/src )

ecbuild_add_test( TARGET  test_nemo_feedback_writer.x
                  SOURCES test_feedback_writer.cc
                  LIBS    nemo_feedback NetCDF::NetCDF_CXX)
set_target_properties(test_nemo_feedback_writer.x PROPERTIES CXX_STANDARD 17)

ecbuild_add_test( TARGET  test_nemo_feedback_data.x
                  SOURCES test_feedback_data.cc
                  LIBS    nemo_feedback)

