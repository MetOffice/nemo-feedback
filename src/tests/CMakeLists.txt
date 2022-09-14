add_subdirectory(testinput)
add_subdirectory(testoutput)
add_subdirectory(nemo-feedback)
add_subdirectory(mains)

ecbuild_add_test( TARGET test_orcamodel_hofx_ice_writer
                  OMP 1
                  ARGS testinput/seaice_obs_writer.yaml
                  COMMAND orcamodel_hofx.x )

ecbuild_add_test( TARGET test_orcamodel_hofx_ice_file
                  OMP 1
                  ARGS testoutput/test_hofx_feedback_writer_out.nc testoutput/test_hofx_feedback_writer_out_ref.cdl
                  COMMAND ./mains/compare_nc_cdl.sh
                  TEST_DEPENDS test_orcamodel_hofx_ice_writer)

ecbuild_add_test( TARGET test_orcamodel_hofx_two_vars_writer
                  OMP 1
                  ARGS testinput/hofx_two_vars_writer.yaml
                  COMMAND orcamodel_hofx.x )

ecbuild_add_test( TARGET test_orcamodel_hofx_two_vars_file
                  OMP 1
                  ARGS testoutput/test_hofx_two_vars_writer_out.nc testoutput/test_hofx_two_vars_writer_out_ref.cdl
                  COMMAND ./mains/compare_nc_cdl.sh
                  TEST_DEPENDS test_orcamodel_hofx_two_vars_writer)

ecbuild_add_test( TARGET test_orcamodel_hofx_profiles_writer
                  OMP 1
                  ARGS testinput/hofx_profiles_writer.yaml
                  COMMAND orcamodel_hofx3D.x )

ecbuild_add_test( TARGET test_orcamodel_hofx_profiles_file
                  OMP 1
                  ARGS testoutput/test_hofx_profiles_writer_out.nc testoutput/test_hofx_profiles_writer_out_ref.cdl
                  COMMAND ./mains/compare_nc_cdl.sh
                  TEST_DEPENDS test_orcamodel_hofx_profiles_writer)

ecbuild_add_test( TARGET test_orcamodel_hofx_profiles_two_vars_writer
                  OMP 1
                  ARGS testinput/hofx_profiles_2var_writer.yaml
                  COMMAND orcamodel_hofx3D.x )

ecbuild_add_test( TARGET test_orcamodel_hofx_profiles_two_vars_file
                  OMP 1
                  ARGS testoutput/test_hofx3d_nc_prof_2vars_writer_out.nc testoutput/test_hofx3d_nc_prof_2vars_writer_out_ref.cdl
                  COMMAND ./mains/compare_nc_cdl.sh
                  TEST_DEPENDS test_orcamodel_hofx_profiles_two_vars_writer)

