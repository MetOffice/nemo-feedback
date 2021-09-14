
list( APPEND fdbk_writer_test_input
  seaice_obs_writer.yaml
)

foreach(FILENAME ${fdbk_writer_test_input})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)
