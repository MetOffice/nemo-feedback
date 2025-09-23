# (C) British Crown Copyright 2025 Met Office

list( APPEND fdbk_writer_test_output
  test_hofx_sic_writer_out_ref.cdl
  test_hofx_two_vars_writer_out_ref.cdl
  test_hofx_profiles_writer_out_ref.cdl
  test_hofx3d_nc_prof_2vars_writer_out_ref.cdl
  test_hofx_mpi_2_writer_out_00001_ref.cdl
  test_hofx_mpi_3_writer_out_00001_ref.cdl
)

foreach(FILENAME ${fdbk_writer_test_output})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)
