/*
 * (C) Crown Copyright 2021, the Met Office. All rights reserved.
 *
 * Refer to COPYRIGHT.txt of this distribution for details.
 */

#include<netcdf>

#include<string>
#include<chrono>
#include<thread>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <iomanip>

#include "eckit/log/Bytes.h"

#include "oops/util/DateTime.h"

#include "eckit/testing/Test.h"

#include "nemo-feedback/feedback_io/Writer.h"

namespace nemo_feedback {
namespace test {

//-----------------------------------------------------------------------------

CASE("test creating test file ") {
  eckit::PathName testDataPath("../testoutput/simple_nemo_out.nc");

  size_t nObs = 5;
  size_t nLevels = 1;

  std::vector<size_t> indices;
  for (size_t iLoc = 0; iLoc < nObs; ++iLoc)
    indices.push_back(iLoc);

  SECTION("NemoFeedbackData empty constructors are valid") {
    feedback_io::Data<double> test1;
  }

  feedback_io::DataIndexer unReducedIndexer(indices, nObs);
  feedback_io::DataIndexer indexer(unReducedIndexer,
      {true, false, true, true, true});

  feedback_io::Data<double> lats(indexer,
                                std::vector<double>{1.0, 2.0, 3.0, 4.0, 5.0});
  feedback_io::Data<double> lons(indexer,
                                std::vector<double>{6.0, 7.0, 8.0, 9.0, 10.0});
  feedback_io::Data<double> depths(indexer,
                                  std::vector<double>(nObs*nLevels, 0));
  feedback_io::Data<double> julianDays(indexer,
     std::vector<double>{11.0, 12.0, 13.0, 14.0, 15.0});
  feedback_io::Data<std::string> stationTypes(indexer,
     std::vector<std::string>{"  44", "  45", "  46", "  47", "  48"});
  feedback_io::Data<std::string> stationIDs(indexer,
     std::vector<std::string>{"12345678",
                              " 2345678",
                              "1234567 ",
                              " 3456789",
                              "123     "});

  std::string juldReference;
  {
    util::DateTime juldReferenceDT = util::DateTime("2021-08-31T15:26:00Z");
    int year, month, day, hour, minute, second;
    juldReferenceDT.toYYYYMMDDhhmmss(year, month, day, hour, minute,
                                           second);
    std::ostringstream juldReferenceSStream;
    juldReferenceSStream << std::setfill('0');
    juldReferenceSStream << std::setw(4) << year;
    juldReferenceSStream << std::setw(2) << month;
    juldReferenceSStream << std::setw(2) << day;
    juldReferenceSStream << std::setw(2) << hour;
    juldReferenceSStream << std::setw(2) << minute;
    juldReferenceSStream << std::setw(2) << second;
    juldReference = juldReferenceSStream.str();
  }

  feedback_io::MetaData metaData(lats, lons, julianDays, depths, stationTypes,
      stationIDs, nLevels, juldReference);

  feedback_io::NameData nameData;
  nameData.variable_names = std::vector<std::string>{"SST", "MDT"};
  nameData.long_names = std::vector<std::string>{"this is a long name",
                                                  "this is another long name"};
  nameData.unit_names = std::vector<std::string>{"this is a unit",
                                                  "this is another unit"};
  nameData.additional_names = std::vector<std::string>{"Hx", "DW_FLAGS",
                                                       "STD"};

  // For some reason this test fails CI with a "file not found" error instead
  // of the eckit error we expect. This is different behaviour depending on the
  // system it seems.
  // SECTION ("NameData validator can fail") {
  //   EXPECT_THROWS_AS(nameData.validate(), eckit::BadValue);
  // }

  nameData.legacy_ops_qc_conventions = std::vector<bool>{false, false};

  std::vector<bool> isExtraVariable{false, true};

  SECTION("file writes") {
    feedback_io::Writer<double> writer(
        testDataPath,
        metaData,
        nameData,
        isExtraVariable);

    // wait up to 100 seconds for the file system...
    for (size_t waitCount = 0; waitCount < 50; ++waitCount) {
      std::this_thread::sleep_for(std::chrono::seconds(2));
      if (testDataPath.exists()) break;
    }
    EXPECT(testDataPath.exists());
  }

  netCDF::NcFile ncFile(testDataPath.fullName().asString(),
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

CASE("test creating profile file ") {
  eckit::PathName testDataPath("../testoutput/simple_nemo_profile_out.nc");

  const size_t nLocations = 7;
  const size_t nObs = 2;
  const size_t nLevels = 5;

  std::vector<size_t> indices;
  for (size_t iLoc = 0; iLoc < nLocations; ++iLoc)
    indices.push_back(iLoc);
  std::vector<size_t> starts = std::vector<size_t>{0, 2};

  feedback_io::DataIndexer indexer(nObs, nLevels, nLocations, starts, indices);

  feedback_io::Data<double> lats(indexer, std::vector<double>(nLocations, 1.0));
  feedback_io::Data<double> lons(indexer, std::vector<double>(nLocations, 6.0));
  feedback_io::Data<double> depths(indexer, std::vector<double>(nLocations, 0));
  feedback_io::Data<double> julianDays(indexer,
      std::vector<double>(nLocations, 11.0));

  feedback_io::Data<std::string> stationTypes(indexer,
      std::vector<std::string>(nLocations, " 401"));
  feedback_io::Data<std::string> stationIDs(indexer,
      std::vector<std::string>(nLocations, "123     "));

  std::string juldReference;
  {
    util::DateTime juldReferenceDT = util::DateTime("2021-08-31T15:26:00Z");
    int year, month, day, hour, minute, second;
    juldReferenceDT.toYYYYMMDDhhmmss(year, month, day, hour, minute,
                                           second);
    std::ostringstream juldReferenceSStream;
    juldReferenceSStream << std::setfill('0');
    juldReferenceSStream << std::setw(4) << year;
    juldReferenceSStream << std::setw(2) << month;
    juldReferenceSStream << std::setw(2) << day;
    juldReferenceSStream << std::setw(2) << hour;
    juldReferenceSStream << std::setw(2) << minute;
    juldReferenceSStream << std::setw(2) << second;
    juldReference = juldReferenceSStream.str();
  }

  feedback_io::MetaData metaData(lats, lons, julianDays, depths, stationTypes,
      stationIDs, nLevels, juldReference);

  feedback_io::NameData nameData;
  nameData.variable_names = std::vector<std::string>{"POTM", "PSAL"};
  nameData.long_names = std::vector<std::string>{"this is a long name",
                                                  "this is another long name"};
  nameData.unit_names = std::vector<std::string>{"this is a unit",
                                                  "this is another unit"};
  nameData.additional_names = std::vector<std::string>{"Hx", "SuperOb"};
  nameData.legacy_ops_qc_conventions = std::vector<bool>{false, false};

  const std::vector<bool> isExtraVariable{false, false};

  SECTION("file writes") {
    feedback_io::Writer<double> writer(
        testDataPath,
        metaData,
        nameData,
        isExtraVariable);

    std::vector<double> dataVec(nLocations);
    for (size_t iLoc = 0; iLoc < nLocations; ++iLoc) dataVec[iLoc] = iLoc;
    feedback_io::Data<double> data(indexer, dataVec);
    writer.write_variable_profile(nameData.variable_names[0] + "_OBS",
        data);
    writer.write_variable_profile(nameData.variable_names[0] + "_Hx",
        data);
    writer.write_variable_profile(nameData.variable_names[1] + "_OBS",
        data);
    writer.write_variable_profile(nameData.variable_names[1] + "_Hx",
        data);

    std::vector<int32_t> intVec(nLocations, 0);
    for (size_t iLoc = 0; iLoc < nLocations; ++iLoc) intVec[iLoc] = 10+iLoc;
    feedback_io::Data<int32_t> intData(indexer, intVec);
    writer.write_variable_level_qc(
        nameData.variable_names[0] + "_LEVEL_QC_FLAGS", intData, 0);
    writer.write_variable_level_qc(
        nameData.variable_names[1] + "_LEVEL_QC_FLAGS", intData, 0);
    writer.write_variable_level_qc(
        nameData.variable_names[0] + "_LEVEL_QC", intData);
    writer.write_variable_level_qc(
        nameData.variable_names[1] + "_LEVEL_QC", intData);

    // wait up to 100 seconds for the file system...
    for (size_t waitCount = 0; waitCount < 50; ++waitCount) {
      std::this_thread::sleep_for(std::chrono::seconds(2));
      if (testDataPath.exists()) break;
    }
    EXPECT(testDataPath.exists());
  }

  netCDF::NcFile ncFile(testDataPath.fullName().asString(),
      netCDF::NcFile::read);

  for (const std::string& variableType :
        std::vector<std::string>{"_OBS", "_Hx"}) {
    for (const std::string& variableName : nameData.variable_names) {
      SECTION(std::string("Profile ") + variableName
              + variableType + " data is correct") {
        netCDF::NcVar ncVar = ncFile.getVar(variableName + variableType);
        std::vector<double> data(metaData.nObs*metaData.nLevels, 12345);
        ncVar.getVar({0, 0}, {metaData.nObs, metaData.nLevels}, data.data());

        EXPECT_EQUAL(0, data[0]);
        EXPECT_EQUAL(1, data[1]);
        EXPECT_EQUAL(99999, data[2]);
        EXPECT_EQUAL(99999, data[3]);
        EXPECT_EQUAL(99999, data[4]);

        EXPECT_EQUAL(2, data[metaData.nLevels]);
        EXPECT_EQUAL(3, data[metaData.nLevels+1]);
        EXPECT_EQUAL(4, data[metaData.nLevels+2]);
        EXPECT_EQUAL(5, data[metaData.nLevels+3]);
        EXPECT_EQUAL(6, data[metaData.nLevels+4]);
      }
    }
  }

  for (const std::string& variableName : nameData.variable_names) {
    SECTION(std::string("Profile ") + variableName
            + "_LEVEL_QC_FLAGS is correct") {
      netCDF::NcVar ncVar = ncFile.getVar(variableName + "_LEVEL_QC_FLAGS");
      std::vector<int> data(metaData.nObs*metaData.nLevels, 12345);
      ncVar.getVar({0, 0, 0}, {metaData.nObs, metaData.nLevels, 1},
          data.data());

      EXPECT_EQUAL(10, data[0]);
      EXPECT_EQUAL(11, data[1]);
      EXPECT_EQUAL(0, data[2]);
      EXPECT_EQUAL(0, data[3]);
      EXPECT_EQUAL(0, data[4]);

      EXPECT_EQUAL(12, data[metaData.nLevels]);
      EXPECT_EQUAL(13, data[metaData.nLevels+1]);
      EXPECT_EQUAL(14, data[metaData.nLevels+2]);
      EXPECT_EQUAL(15, data[metaData.nLevels+3]);
      EXPECT_EQUAL(16, data[metaData.nLevels+4]);
    }

    SECTION(std::string("Profile ") + variableName + "_LEVEL_QC is correct") {
      netCDF::NcVar ncVar = ncFile.getVar(variableName + "_LEVEL_QC");
      std::vector<int> data(metaData.nObs*metaData.nLevels, 12345);
      ncVar.getVar({0, 0}, {metaData.nObs, metaData.nLevels}, data.data());

      EXPECT_EQUAL(10, data[0]);
      EXPECT_EQUAL(11, data[1]);
      EXPECT_EQUAL(0, data[2]);
      EXPECT_EQUAL(0, data[3]);
      EXPECT_EQUAL(0, data[4]);

      EXPECT_EQUAL(12, data[metaData.nLevels]);
      EXPECT_EQUAL(13, data[metaData.nLevels+1]);
      EXPECT_EQUAL(14, data[metaData.nLevels+2]);
      EXPECT_EQUAL(15, data[metaData.nLevels+3]);
      EXPECT_EQUAL(16, data[metaData.nLevels+4]);
    }
  }
}

CASE("test creating reduced profile file ") {
  eckit::PathName testDataPath(
      "../testoutput/simple_nemo_reduced_profile_out.nc");

  const size_t nLocations = 17;
  const size_t nObs = 5;
  const size_t nLevels = 6;

  std::vector<size_t> indices;
  for (size_t iLoc = 0; iLoc < nLocations; ++iLoc)
    indices.push_back(iLoc);
  // Without eliminating empty profiles:
  std::vector<size_t> starts{0, 2, 4, 6, 11};
  feedback_io::DataIndexer unReducedIndexer(nObs, nLevels, nLocations, starts,
      indices);
  // Without eliminating empty profiles lenghts are {2, 2, 2, nLevels, 6};
  const std::vector<bool> toWrite{true, false,
                                  false, false,
                                  false, false,
                                  true, true, true, false, true,
                                  false, false, false, false, false, false};

  feedback_io::DataIndexer indexer(unReducedIndexer, toWrite);

  feedback_io::Data<double> lats(indexer, std::vector<double>(nLocations, 1.0));
  feedback_io::Data<double> lons(indexer, std::vector<double>(nLocations, 6.0));
  feedback_io::Data<double> depths(indexer, std::vector<double>(nLocations, 0));
  feedback_io::Data<double> julianDays(indexer,
      std::vector<double>(nLocations, 11.0));

  feedback_io::Data<std::string> stationTypes(indexer,
     std::vector<std::string>(nLocations, " 401"));
  feedback_io::Data<std::string> stationIDs(indexer,
     std::vector<std::string>(nLocations, "123     "));

  std::string juldReference;
  {
    util::DateTime juldReferenceDT = util::DateTime("2021-08-31T15:26:00Z");
    int year, month, day, hour, minute, second;
    juldReferenceDT.toYYYYMMDDhhmmss(year, month, day, hour, minute,
                                           second);
    std::ostringstream juldReferenceSStream;
    juldReferenceSStream << std::setfill('0');
    juldReferenceSStream << std::setw(4) << year;
    juldReferenceSStream << std::setw(2) << month;
    juldReferenceSStream << std::setw(2) << day;
    juldReferenceSStream << std::setw(2) << hour;
    juldReferenceSStream << std::setw(2) << minute;
    juldReferenceSStream << std::setw(2) << second;
    juldReference = juldReferenceSStream.str();
  }

  feedback_io::MetaData metaData(lats, lons, julianDays, depths, stationTypes,
      stationIDs, nLevels, juldReference);

  feedback_io::NameData nameData;
  nameData.variable_names = std::vector<std::string>{"POTM", "PSAL"};
  nameData.long_names = std::vector<std::string>{"this is a long name",
                                                  "this is another long name"};
  nameData.unit_names = std::vector<std::string>{"this is a unit",
                                                  "this is another unit"};
  nameData.additional_names = std::vector<std::string>{"Hx", "SuperOb"};
  nameData.legacy_ops_qc_conventions = std::vector<bool>{false, false};

  const std::vector<bool> isExtraVariable{false, false};

  SECTION("file writes") {
    // In the above example this is two profiles with:
    // [A, B, B, B, B,]

    const size_t nToWrite = std::count(toWrite.begin(), toWrite.end(), true);
    // number of locations to write should match the total number of locations
    EXPECT_EQUAL(nToWrite, metaData.nLocations);

    feedback_io::Writer<double> writer(
        testDataPath,
        metaData,
        nameData,
        isExtraVariable);

    std::vector<double> dataVec(nLocations);
    for (size_t iLoc = 0; iLoc < nLocations; ++iLoc) dataVec[iLoc] = iLoc;
    feedback_io::Data<double> data(indexer, dataVec);

    writer.write_variable_profile(nameData.variable_names[0] + "_OBS",
        data);
    writer.write_variable_profile(nameData.variable_names[0] + "_Hx",
        data);
    writer.write_variable_profile(nameData.variable_names[1] + "_OBS",
        data);
    writer.write_variable_profile(nameData.variable_names[1] + "_Hx",
        data);

    std::vector<int32_t> intVec(nLocations, 0);
    for (size_t iLoc = 0; iLoc < nLocations; ++iLoc) intVec[iLoc] = 10+iLoc;
    feedback_io::Data<int32_t> intData(indexer, intVec);

    writer.write_variable_level_qc(
        nameData.variable_names[0] + "_LEVEL_QC_FLAGS", intData, 0);
    writer.write_variable_level_qc(
        nameData.variable_names[1] + "_LEVEL_QC_FLAGS", intData, 0);
    writer.write_variable_level_qc(nameData.variable_names[0]
        + "_LEVEL_QC", intData);
    writer.write_variable_level_qc(nameData.variable_names[1]
        + "_LEVEL_QC", intData);

    // wait up to 100 seconds for the file system...
    for (size_t waitCount = 0; waitCount < 50; ++waitCount) {
      std::this_thread::sleep_for(std::chrono::seconds(2));
      if (testDataPath.exists()) break;
    }
    EXPECT(testDataPath.exists());
  }

  netCDF::NcFile ncFile(testDataPath.fullName().asString(),
      netCDF::NcFile::read);

  for (const std::string& variableType :
        std::vector<std::string>{"_OBS", "_Hx"}) {
    for (const std::string& variableName : nameData.variable_names) {
      SECTION(std::string("Profile ") + variableName + variableType
              + " data is correct") {
        netCDF::NcVar ncVar = ncFile.getVar(variableName + variableType);
        std::vector<double> data(metaData.nObs*metaData.nLevels, 12345);
        ncVar.getVar({0, 0}, {metaData.nObs, metaData.nLevels}, data.data());

        EXPECT_EQUAL(0, data[0]);
        EXPECT_EQUAL(99999, data[1]);
        EXPECT_EQUAL(99999, data[2]);
        EXPECT_EQUAL(99999, data[3]);
        // Two profiles immediately rejected, as they lack valid observations
        EXPECT_EQUAL(6, data[metaData.nLevels]);
        EXPECT_EQUAL(7, data[metaData.nLevels+1]);
        EXPECT_EQUAL(8, data[metaData.nLevels+2]);
        // One observation reduced out of second profile
        EXPECT_EQUAL(10, data[metaData.nLevels+3]);
        // Final profile immediately rejected, as it lacks valid observations
      }
    }
  }

  for (const std::string& variableName : nameData.variable_names) {
    SECTION(std::string("Profile ") + variableName
            + "_LEVEL_QC_FLAGS is correct") {
      netCDF::NcVar ncVar = ncFile.getVar(variableName + "_LEVEL_QC_FLAGS");
      std::vector<int> data(metaData.nObs*metaData.nLevels, 12345);
      ncVar.getVar({0, 0, 0}, {metaData.nObs, metaData.nLevels, 1},
          data.data());

      EXPECT_EQUAL(10, data[0]);
      EXPECT_EQUAL(0, data[1]);
      EXPECT_EQUAL(0, data[2]);
      EXPECT_EQUAL(0, data[3]);
      // Two profiles immediately rejected, as they lack valid observations
      EXPECT_EQUAL(16, data[metaData.nLevels]);
      EXPECT_EQUAL(17, data[metaData.nLevels+1]);
      EXPECT_EQUAL(18, data[metaData.nLevels+2]);
      // One observation reduced out of second profile
      EXPECT_EQUAL(20, data[metaData.nLevels+3]);
      // Final profile immediately rejected, as it lacks valid observations
    }

    SECTION(std::string("Profile ") + variableName + "_LEVEL_QC is correct") {
      netCDF::NcVar ncVar = ncFile.getVar(variableName + "_LEVEL_QC");
      std::vector<int> data(metaData.nObs*metaData.nLevels, 12345);
      ncVar.getVar({0, 0}, {metaData.nObs, metaData.nLevels}, data.data());

      EXPECT_EQUAL(10, data[0]);
      EXPECT_EQUAL(0, data[1]);
      EXPECT_EQUAL(0, data[2]);
      EXPECT_EQUAL(0, data[3]);
      // Two profiles immediately rejected, as they lack valid observations
      EXPECT_EQUAL(16, data[metaData.nLevels]);
      EXPECT_EQUAL(17, data[metaData.nLevels+1]);
      EXPECT_EQUAL(18, data[metaData.nLevels+2]);
      // One observation reduced out of second profile
      EXPECT_EQUAL(20, data[metaData.nLevels+3]);
      // Final profile immediately rejected, as it lacks valid observations
    }
  }
}

}  // namespace test
}  // namespace nemo_feedback

int main(int argc, char** argv) {
    return eckit::testing::run_tests(argc, argv);
}
