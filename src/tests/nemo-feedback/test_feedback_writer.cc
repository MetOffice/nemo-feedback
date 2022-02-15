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

  size_t n_obs = 5;
  size_t n_levels = 1;
  std::vector<double> lats{1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<double> lons{6.0, 7.0, 8.0, 9.0, 10.0};
  std::vector<double> depths(n_obs*n_levels, 0);
  std::vector<double> times{11.0, 12.0, 13.0, 14.0, 15.0};
  std::vector<std::string> variable_names{"SST"};
  std::vector<std::string> long_names{"this is a long name"};
  std::vector<std::string> unit_names{"this is a unit"};
  std::vector<std::string> additional_variables{"Hx", "DW_FLAGS", "STD"};
  util::DateTime juld_reference("2021-08-31T15:26:00Z");
  std::vector<std::string> station_types{"  44", "  45", "  46", "  47", "  48"};

  SECTION("file writes") {
    NemoFeedbackWriter fdbk_writer(
        test_data_path, 
        4,
        {true, false, true, true, true},
        lons, 
        lats, 
        depths,
        times, 
        variable_names, 
        long_names, 
        unit_names,
        additional_variables, 
        n_levels, 
        juld_reference, 
        station_types);
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
  
  SECTION("STATION_TYPE variable is correct") {
    netCDF::NcVar ncVar = ncFile.getVar("STATION_TYPE");
    char data[4] = { ' ' };
    ncVar.getVar({0, 0}, {1, 4}, data);
    EXPECT_EQUAL(static_cast<std::string>(data).substr(0, 4),
        static_cast<std::string>("  44"));
    ncVar.getVar({1, 0}, {1, 4}, data);
    EXPECT_EQUAL(static_cast<std::string>(data).substr(0, 4),
        static_cast<std::string>("  46"));
    ncVar.getVar({2, 0}, {1, 4}, data);
    EXPECT_EQUAL(static_cast<std::string>(data).substr(0, 4),
        static_cast<std::string>("  47"));
    ncVar.getVar({3, 0}, {1, 4}, data);
    EXPECT_EQUAL(static_cast<std::string>(data).substr(0, 4),
        static_cast<std::string>("  48"));
  }
}

}  // namespace test
}  // namespace nemo_feedback

int main(int argc, char** argv) {
    return eckit::testing::run_tests(argc, argv);
}
