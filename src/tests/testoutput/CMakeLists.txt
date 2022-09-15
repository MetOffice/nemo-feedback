
list( APPEND fdbk_writer_test_output
  test_hofx_feedback_writer.ref
  test_hofx_feedback_writer_out_ref.cdl
  test_hofx_two_vars_writer.ref
  test_hofx_two_vars_writer_out_ref.cdl
  test_hofx_profiles_writer.ref
  test_hofx_profiles_writer_out_ref.cdl
  test_hofx3d_nc_prof_2vars_writer.ref
  test_hofx3d_nc_prof_2vars_writer_out_ref.cdl
)

foreach(FILENAME ${fdbk_writer_test_output})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)
