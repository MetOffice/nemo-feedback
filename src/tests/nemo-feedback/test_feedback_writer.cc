/*
 * (C) Crown Copyright 2021, the Met Office. All rights reserved.
 *
 * Refer to COPYRIGHT.txt of this distribution for details. 
 */

#include "eckit/log/Bytes.h"

#include "oops/util/DateTime.h"

#include "eckit/testing/Test.h"

#include "nemo-feedback/NemoFeedbackWriter.h"

namespace nemo_feedback {
namespace test {

//-----------------------------------------------------------------------------

CASE("test creating test file ") {
  eckit::PathName test_data_path("simple_nemo_out.nc");

  size_t n_obs = 3;
  size_t n_levels = 1;
  std::vector<double> lats(n_obs, 0);
  std::vector<double> lons(n_obs, 0);
  std::vector<double> depths(n_obs*n_levels, 0);
  std::vector<double> times(n_obs, 0);
  std::vector<std::string> variable_names{"SST"};
  std::vector<std::string> additional_variables{"Hx", "DW_FLAGS", "STD"};
  util::DateTime juld_reference("2021-08-31T15:26:00Z");

  NemoFeedbackWriter fdbk_writer(test_data_path, lons, lats, depths,
      times, variable_names, additional_variables, n_levels, juld_reference);
}

}  // namespace test
}  // namespace nemo_feedback

int main(int argc, char** argv) {
    return eckit::testing::run_tests(argc, argv);
}
