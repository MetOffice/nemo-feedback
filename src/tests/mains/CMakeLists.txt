ecbuild_add_executable( TARGET  test_NemoInputsFilters.x
                        SOURCES TestObsFilters.cc
                        LIBS    nemo_feedback ufo oops)

list( APPEND script_files
  compare_nc_cdl.sh
)

foreach(FILENAME ${script_files})
execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
                ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)

