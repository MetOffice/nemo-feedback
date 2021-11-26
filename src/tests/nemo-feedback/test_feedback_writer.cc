/*
 * (C) Crown Copyright 2021, the Met Office. All rights reserved.
 *
 * Refer to COPYRIGHT.txt of this distribution for details.
 */

#include<netcdf>
#include<string>

#include "eckit/log/Bytes.h"

#include "oops/util/DateTime.h"

#include "eckit/testing/Test.h"

#include "nemo-feedback/NemoFeedbackWriter.h"

namespace nemo_feedback {
namespace test {

//-----------------------------------------------------------------------------

CASE("test creating test file ") {
  eckit::PathName test_data_path("../testoutput/simple_nemo_out.nc");

  size_t n_obs = 3;
  size_t n_levels = 1;
  std::vector<double> lats(n_obs, 0);
  std::vector<double> lons(n_obs, 0);
  std::vector<double> depths(n_obs*n_levels, 0);
  std::vector<double> times(n_obs, 0);
  std::vector<std::string> variable_names{"SST"};
  std::vector<std::string> additional_variables{"Hx", "DW_FLAGS", "STD"};
  util::DateTime juld_reference("2021-08-31T15:26:00Z");

  SECTION("file writes") {
    NemoFeedbackWriter fdbk_writer(test_data_path, lons, lats, depths,
        times, variable_names, additional_variables, n_levels, juld_reference);
  }

  netCDF::NcFile ncFile(test_data_path.fullName().asString(),
      netCDF::NcFile::read);

  SECTION("JULD_REFERENCE variable is correct") {
    netCDF::NcVar ncVar = ncFile.getVar("JULD_REFERENCE");
    char data[15] = { ' ' };
    ncVar.getVar({0}, {14}, data);

    EXPECT_EQUAL(static_cast<std::string>(data).substr(0, 14),
        static_cast<std::string>("20210831152600"));
  }

  SECTION("ENTRIES variable is correct") {
    netCDF::NcVar entriesVar = ncFile.getVar("ENTRIES");
    char data[9] = { ' ' };
    entriesVar.getVar({0, 0}, {1, 8}, data);
    EXPECT_EQUAL(static_cast<std::string>(data).substr(0, 8),
        static_cast<std::string>("Hx      "));
    entriesVar.getVar({1, 0}, {1, 8}, data);
    EXPECT_EQUAL(static_cast<std::string>(data).substr(0, 8),
        static_cast<std::string>("DW_FLAGS"));
    entriesVar.getVar({2, 0}, {1, 8}, data);
    EXPECT_EQUAL(static_cast<std::string>(data).substr(0, 8),
        static_cast<std::string>("STD     "));
  }

  SECTION("VARIABLES variable is correct") {
    netCDF::NcVar ncVar = ncFile.getVar("VARIABLES");
    char data[9] = { ' ' };
    ncVar.getVar({0, 0}, {1, 8}, data);

    EXPECT_EQUAL(static_cast<std::string>(data).substr(0, 8),
        static_cast<std::string>("SST     "));
  }
}

}  // namespace test
}  // namespace nemo_feedback

int main(int argc, char** argv) {
    return eckit::testing::run_tests(argc, argv);
}
