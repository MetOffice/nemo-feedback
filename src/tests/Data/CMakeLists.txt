# (C) British Crown Copyright 2024 Met Office

list( APPEND nemofeedback_test_data
hofx_sic_obs.nc
hofx_two_vars_obs.nc
hofx_potm_obs.nc
hofx_prof_2var_obs.nc
test_hofx3d_nc_prof_2vars_writer_geovals.nc
test_hofx_profiles_writer_geovals.nc
test_hofx_sic_writer_geovals.nc
test_hofx_two_vars_writer_geovals.nc
)

foreach(FILENAME ${nemofeedback_test_data})
    execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
           ${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}
           ${CMAKE_CURRENT_BINARY_DIR}/${FILENAME} )
endforeach(FILENAME)
