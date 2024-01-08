# (C) British Crown Copyright 2024 Met Office

list( APPEND fdbk_writer_test_input
  hofx_sic_writer.yaml
  hofx_two_vars_writer.yaml
  hofx_profiles_writer.yaml
  hofx_profiles_2var_writer.yaml
  test_name_map.yaml
)

foreach(FILENAME ${fdbk_writer_test_input})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)
