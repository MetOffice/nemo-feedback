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

  NemoFeedbackWriter::CoordData coords;
  coords.n_levels = 1;
  coords.lats = std::vector<double>{1.0, 2.0, 3.0, 4.0, 5.0};
  coords.lons = std::vector<double>{6.0, 7.0, 8.0, 9.0, 10.0};
  coords.depths = std::vector<double>(n_obs*coords.n_levels, 0);
  coords.julian_days = std::vector<double>{11.0, 12.0, 13.0, 14.0, 15.0};
  coords.juld_reference = util::DateTime("2021-08-31T15:26:00Z");

  NemoFeedbackWriter::NameData name_data;
  name_data.variable_names = std::vector<std::string>{"SST", "MDT"};
  name_data.long_names = std::vector<std::string>{"this is a long name",
                                                  "this is another long name"};
  name_data.unit_names = std::vector<std::string>{"this is a unit",
                                                  "this is another unit"};
  name_data.additional_names = std::vector<std::string>{"Hx", "DW_FLAGS",
                                                       "STD"};

  std::vector<bool> extra_variables{false, true};
  std::vector<std::string> station_types{"  44", "  45", "  46", "  47",
                                         "  48"};
  std::vector<std::string> station_ids{"12345678",
                                       " 2345678",
                                       "1234567 ",
                                       " 3456789",
                                       "123     "};

  SECTION("file writes") {
    NemoFeedbackWriter fdbk_writer(
        test_data_path,
        4,
        {true, false, true, true, true},
        coords,
        name_data,
        extra_variables,
        station_types,
        station_ids);

    EXPECT(test_data_path.exists());
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
    char data[5] = { ' ' };
    data[4] = '\0';
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

  SECTION("STATION_IDENTIFIER variable is correct") {
    netCDF::NcVar ncVar = ncFile.getVar("STATION_IDENTIFIER");
    char data[9] = { ' ' };
    data[8] = '\0';
    ncVar.getVar({0, 0}, {1, 8}, data);
    EXPECT_EQUAL(static_cast<std::string>(data).substr(0, 8),
        static_cast<std::string>("12345678"));
    ncVar.getVar({1, 0}, {1, 8}, data);
    EXPECT_EQUAL(static_cast<std::string>(data).substr(0, 8),
        static_cast<std::string>("1234567 "));
    ncVar.getVar({2, 0}, {1, 8}, data);
    EXPECT_EQUAL(static_cast<std::string>(data).substr(0, 8),
        static_cast<std::string>(" 3456789"));
    ncVar.getVar({3, 0}, {1, 8}, data);
    EXPECT_EQUAL(static_cast<std::string>(data).substr(0, 8),
        static_cast<std::string>("123     "));
  }
}

CASE("test creating test profile file ") {
  eckit::PathName test_data_path("../testoutput/simple_nemo_profile_out.nc");

  const size_t n_locs = 7;
  const size_t n_obs = 2;
  NemoFeedbackWriter::CoordData coords;
  coords.n_levels = 5;
  coords.lats = std::vector<double>(n_obs, 1.0);
  coords.lons = std::vector<double>(n_obs, 6.0);
  coords.depths = std::vector<double>(n_locs, 0);
  coords.julian_days = std::vector<double>(n_obs, 11.0);
  coords.juld_reference = util::DateTime("2021-08-31T15:26:00Z");

  NemoFeedbackWriter::NameData name_data;
  name_data.variable_names = std::vector<std::string>{"POTM", "PSAL"};
  name_data.long_names = std::vector<std::string>{"this is a long name",
                                                  "this is another long name"};
  name_data.unit_names = std::vector<std::string>{"this is a unit",
                                                  "this is another unit"};
  name_data.additional_names = std::vector<std::string>{"Hx", "SuperOb"};

  const std::vector<bool> extra_variables{false, false};
  const std::vector<std::string> station_types{" 401", " 401"};
  const std::vector<std::string> station_ids{"123     ", "123     "};

  SECTION("file writes") {
    NemoFeedbackWriter fdbk_writer(
        test_data_path,
        n_obs,
        {true, true},
        coords,
        name_data,
        extra_variables,
        station_types,
        station_ids);

    const std::vector<size_t> record_starts{0, 2}, record_counts{2, coords.n_levels};

    std::vector<double> data(n_locs, 0);
    for (int i = 0; i < n_locs; ++i) data[i] = i;
    fdbk_writer.write_variable_profile(name_data.variable_names[0] + "_OBS",
        data, record_starts, record_counts);
    fdbk_writer.write_variable_profile(name_data.variable_names[0] + "_Hx",
        data, record_starts, record_counts);
    fdbk_writer.write_variable_profile(name_data.variable_names[1] + "_OBS",
        data, record_starts, record_counts);
    fdbk_writer.write_variable_profile(name_data.variable_names[1] + "_Hx",
        data, record_starts, record_counts);

    std::vector<int32_t> int_data(n_locs, 0);
    for (int i = 0; i < n_locs; ++i) int_data[i] = 10+i;
    fdbk_writer.write_variable_level_qc(
        name_data.variable_names[0] + "_LEVEL_QC_FLAGS", int_data, 0,
        record_starts, record_counts);
    fdbk_writer.write_variable_level_qc(
        name_data.variable_names[1] + "_LEVEL_QC_FLAGS", int_data, 0,
        record_starts, record_counts);
    fdbk_writer.write_variable_level_qc(
        name_data.variable_names[0] + "_LEVEL_QC", int_data,
        record_starts, record_counts);
    fdbk_writer.write_variable_level_qc(
        name_data.variable_names[1] + "_LEVEL_QC", int_data,
        record_starts, record_counts);

    EXPECT(test_data_path.exists());
  }

  netCDF::NcFile ncFile(test_data_path.fullName().asString(),
      netCDF::NcFile::read);

  for (const std::string& v_type : std::vector<std::string>{"_OBS", "_Hx"}) {
    for (const std::string& v_name : name_data.variable_names) {
      SECTION(std::string("Profile ") + v_name + v_type + " data is correct") {
        netCDF::NcVar ncVar = ncFile.getVar(v_name + v_type);
        std::vector<double> data(n_obs*coords.n_levels, 12345);
        ncVar.getVar({0, 0}, {n_obs, coords.n_levels}, data.data());

        EXPECT_EQUAL(0, data[0]);
        EXPECT_EQUAL(1, data[1]);
        EXPECT_EQUAL(99999, data[2]);
        EXPECT_EQUAL(99999, data[3]);
        EXPECT_EQUAL(99999, data[4]);

        EXPECT_EQUAL(2, data[coords.n_levels]);
        EXPECT_EQUAL(3, data[coords.n_levels+1]);
        EXPECT_EQUAL(4, data[coords.n_levels+2]);
        EXPECT_EQUAL(5, data[coords.n_levels+3]);
        EXPECT_EQUAL(6, data[coords.n_levels+4]);
      }
    }
  }

  for (const std::string& v_name : name_data.variable_names) {
    SECTION(std::string("Profile ") + v_name + "_LEVEL_QC_FLAGS is correct") {
      netCDF::NcVar ncVar = ncFile.getVar(v_name + "_LEVEL_QC_FLAGS");
      std::vector<int> data(n_obs*coords.n_levels, 12345);
      ncVar.getVar({0, 0, 0}, {n_obs, coords.n_levels, 1}, data.data());

      EXPECT_EQUAL(10, data[0]);
      EXPECT_EQUAL(11, data[1]);
      EXPECT_EQUAL(0, data[2]);
      EXPECT_EQUAL(0, data[3]);
      EXPECT_EQUAL(0, data[4]);

      EXPECT_EQUAL(12, data[coords.n_levels]);
      EXPECT_EQUAL(13, data[coords.n_levels+1]);
      EXPECT_EQUAL(14, data[coords.n_levels+2]);
      EXPECT_EQUAL(15, data[coords.n_levels+3]);
      EXPECT_EQUAL(16, data[coords.n_levels+4]);
    }

    SECTION(std::string("Profile ") + v_name + "_LEVEL_QC is correct") {
      netCDF::NcVar ncVar = ncFile.getVar(v_name + "_LEVEL_QC");
      std::vector<int> data(n_obs*coords.n_levels, 12345);
      ncVar.getVar({0, 0}, {n_obs, coords.n_levels}, data.data());

      EXPECT_EQUAL(10, data[0]);
      EXPECT_EQUAL(11, data[1]);
      EXPECT_EQUAL(0, data[2]);
      EXPECT_EQUAL(0, data[3]);
      EXPECT_EQUAL(0, data[4]);

      EXPECT_EQUAL(12, data[coords.n_levels]);
      EXPECT_EQUAL(13, data[coords.n_levels+1]);
      EXPECT_EQUAL(14, data[coords.n_levels+2]);
      EXPECT_EQUAL(15, data[coords.n_levels+3]);
      EXPECT_EQUAL(16, data[coords.n_levels+4]);
    }
  }

}

}  // namespace test
}  // namespace nemo_feedback

int main(int argc, char** argv) {
    return eckit::testing::run_tests(argc, argv);
}
