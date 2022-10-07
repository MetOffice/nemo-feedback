
list( APPEND orcajedi_test_data
hofx_sic_obs.nc
hofx_two_vars_obs.nc
hofx_potm_obs.nc
hofx_prof_2var_obs.nc
)

foreach(FILENAME ${orcajedi_test_data})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)
