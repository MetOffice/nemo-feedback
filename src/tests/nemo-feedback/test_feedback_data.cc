/*
 * (C) British Crown Copyright 2024 Met Office
 */

#include <iostream>
#include <memory>
#include <string>
#include <algorithm>
#include <numeric>

#include "eckit/log/Bytes.h"

#include "oops/util/DateTime.h"

#include "eckit/testing/Test.h"

#include "nemo-feedback/feedback_io/DataIndexer.h"
#include "nemo-feedback/feedback_io/Data.h"

namespace nemo_feedback {
namespace feedback_io {
namespace test {

//-----------------------------------------------------------------------------

/// \brief test cases fixture
static struct Fixture {
  Fixture() {
    // test cases for the vector of indices
    indexingTypes.push_back   ("ordered");  // NOLINT(*)
    indices.emplace_back  (std::vector<size_t>{0, 1, 2, 3, 4, 5, 6, 7, 8});  // NOLINT(*)
    unReducedData.emplace_back(std::vector<size_t>{1, 2, 3, 4, 5, 6, 7, 8, 9});  // NOLINT(*)
    indexingTypes.push_back   ("reversed");  // NOLINT(*)
    indices.emplace_back  (std::vector<size_t>{8, 7, 6, 5, 4, 3, 2, 1, 0});  // NOLINT(*)
    unReducedData.emplace_back(std::vector<size_t>{9, 8, 7, 6, 5, 4, 3, 2, 1});  // NOLINT(*)
    indexingTypes.push_back   ("alternating");  // NOLINT(*)
    indices.emplace_back  (std::vector<size_t>{8, 1, 6, 3, 4, 5, 2, 7, 0});  // NOLINT(*)
    unReducedData.emplace_back(std::vector<size_t>{9, 2, 7, 4, 5, 6, 3, 8, 1});  // NOLINT(*)

    ReducedData.resize(indexingTypes.size());
    nProfiles.resize(indexingTypes.size());
    nLevels.resize(indexingTypes.size());

    // test cases for the selection of data
    selectionTypes.push_back       ("all");  // NOLINT(*)
    selections.emplace_back    (std::vector<bool  >{1, 1, 1, 1, 1, 1, 1, 1, 1});  // NOLINT(*)
    ReducedData[0].emplace_back(std::vector<size_t>{1, 2, 3, 4, 5, 6, 7, 8, 9});  // NOLINT(*)
    nProfiles[0].push_back     (3);  // NOLINT(*)
    nLevels[0].push_back       (4);  // NOLINT(*)
    ReducedData[1].emplace_back(std::vector<size_t>{9, 8, 7, 6, 5, 4, 3, 2, 1});  // NOLINT(*)
    nProfiles[1].push_back     (3);  // NOLINT(*)
    nLevels[1].push_back       (4);  // NOLINT(*)
    ReducedData[2].emplace_back(std::vector<size_t>{9, 2, 7, 4, 5, 6, 3, 8, 1});  // NOLINT(*)
    nProfiles[2].push_back     (3);  // NOLINT(*)
    nLevels[2].push_back       (4);  // NOLINT(*)

    selectionTypes.push_back       ("none");  // NOLINT(*)
    selections.emplace_back    (std::vector<bool  >{0, 0, 0, 0, 0, 0, 0, 0, 0});  // NOLINT(*)
    ReducedData[0].emplace_back(std::vector<size_t>{                         });  // NOLINT(*)
    nProfiles[0].push_back     (0);  // NOLINT(*)
    nLevels[0].push_back       (0);  // NOLINT(*)
    ReducedData[1].emplace_back(std::vector<size_t>{                         });  // NOLINT(*)
    nProfiles[1].push_back     (0);  // NOLINT(*)
    nLevels[1].push_back       (0);  // NOLINT(*)
    ReducedData[2].emplace_back(std::vector<size_t>{                         });  // NOLINT(*)
    nProfiles[2].push_back     (0);  // NOLINT(*)
    nLevels[2].push_back       (0);  // NOLINT(*)

    selectionTypes.push_back       ("subset 1");  // NOLINT(*)
    selections.emplace_back    (std::vector<bool  >{0, 1, 1, 1, 0, 0, 1, 0, 1});  // NOLINT(*)
    ReducedData[0].emplace_back(std::vector<size_t>{   2, 3, 4,       7,    9});  // NOLINT(*)
    nProfiles[0].push_back     (3);  // NOLINT(*)
    nLevels[0].push_back       (2);  // NOLINT(*)
    ReducedData[1].emplace_back(std::vector<size_t>{9,    7,       4, 3, 2,  });  // NOLINT(*)
    nProfiles[1].push_back     (3);  // NOLINT(*)
    nLevels[1].push_back       (2);  // NOLINT(*)
    ReducedData[2].emplace_back(std::vector<size_t>{9, 2, 7, 4,       3,     });  // NOLINT(*)
    nProfiles[2].push_back     (3);  // NOLINT(*)
    nLevels[2].push_back       (2);  // NOLINT(*)

    selectionTypes.push_back       ("subset 2");  // NOLINT(*)
    selections.emplace_back    (std::vector<bool  >{0, 1, 1, 1, 0, 1, 0, 0, 0});  // NOLINT(*)
    ReducedData[0].emplace_back(std::vector<size_t>{   2, 3, 4,    6,        });  // NOLINT(*)
    nProfiles[0].push_back     (2);  // NOLINT(*)
    nLevels[0].push_back       (3);  // NOLINT(*)
    ReducedData[1].emplace_back(std::vector<size_t>{         6,    4, 3, 2   });  // NOLINT(*)
    nProfiles[1].push_back     (2);  // NOLINT(*)
    nLevels[1].push_back       (2);  // NOLINT(*)
    ReducedData[2].emplace_back(std::vector<size_t>{   2,    4,    6, 3,     });  // NOLINT(*)
    nProfiles[2].push_back     (3);  // NOLINT(*)
    nLevels[2].push_back       (2);  // NOLINT(*)
  }

  const std::vector<size_t>      data       = {1, 2, 3, 4, 5, 6, 7, 8, 9};  // NOLINT(*)
  const std::vector<float>       floatData  = {1, 2, 3, 4, 5, 6, 7, 8, 9};  // NOLINT(*)
  const std::vector<double>      doubleData = {1, 2, 3, 4, 5, 6, 7, 8, 9};  // NOLINT(*)
  const std::vector<int32_t>     intData    = {1, 2, 3, 4, 5, 6, 7, 8, 9};  // NOLINT(*)
  const std::vector<std::string> stringData = {"A", "B", "C", "D", "E",  // NOLINT(*)
                                               "F", "G", "H", "I"};
  // imagine our test data is also split into 3 profiles, with 4 max levels
  const size_t unreducedNProfiles = 3;
  const size_t unreducedNLevels = 4;
  const std::vector<size_t> starts = {0, 2, 6};
  std::vector<std::string> indexingTypes;
  std::vector<std::string> selectionTypes;
  std::vector<std::vector<size_t>> indices;
  std::vector<std::vector<bool>> selections;
  std::vector<std::vector<size_t>> unReducedData;

  std::vector<std::vector<std::vector<size_t>>> ReducedData;
  std::vector<std::vector<size_t>> nProfiles;
  std::vector<std::vector<size_t>> nLevels;
} fixture;

CASE("Test Nemo Feedback indexer") {
  for (size_t test = 0; test < fixture.indices.size(); ++test) {
    std::string indexingType = fixture.indexingTypes.at(test);
    std::vector<size_t> indices = fixture.indices.at(test);

    DataIndexer indexer(indices, fixture.unReducedData[test].size());
    SECTION("indexing surface " + indexingType + " indices") {
      for (size_t iOb = 0; iOb < indexer.n_obs(); ++iOb) {
        EXPECT_EQUAL(fixture.unReducedData.at(test).at(iOb),
            fixture.data.at(indexer.at(iOb)));
      }
    }

    for (size_t iSelect = 0; iSelect < fixture.selections.size();
        ++iSelect) {
      SECTION("reduced indexing of surface " + indexingType + " indices with "
              + fixture.selectionTypes[iSelect] + " selected") {
        std::vector<bool> select = fixture.selections[iSelect];
        std::vector<size_t> reducedData =
          fixture.ReducedData[test][iSelect];

        DataIndexer reducedIndexer(indexer, select);
        EXPECT_EQUAL(reducedIndexer.n_obs(), reducedData.size());
        std::cout << "reducedIndexer indices: ";
        for (size_t iOb = 0; iOb < reducedIndexer.n_obs(); ++iOb)
          std::cout << iOb << ": " << reducedIndexer(iOb) << " ";
        std::cout << std::endl;
        size_t iObReduced = 0;
        for (size_t iOb = 0; iOb < indexer.n_obs(); ++iOb) {
          if (select[indexer(iOb)]) {
            EXPECT_EQUAL(reducedData[iObReduced],
                fixture.data[reducedIndexer(iObReduced)]);
            iObReduced++;
          }
        }
      }
    }
  }
  for (size_t test = 0; test < fixture.indices.size(); ++test) {
    std::string indexingType = fixture.indexingTypes[test];
    std::vector<size_t> indices = fixture.indices[test];

    DataIndexer indexer(fixture.unreducedNProfiles,
        fixture.unreducedNLevels, fixture.unReducedData[test].size(),
        fixture.starts, indices);

    SECTION("indexing profile " + indexingType + " indices") {
      size_t iLoc = 0;
      for (size_t iProf = 0; iProf < indexer.n_obs(); ++iProf) {
        for (size_t iLevel = 0; iLevel < indexer.length(iProf); ++iLevel) {
          EXPECT_EQUAL(fixture.unReducedData[test][iLoc],
                       fixture.data[indexer.at(iProf, iLevel)]);
          iLoc++;
        }
      }
      for (size_t iProf = 0; iProf < indexer.n_obs(); ++iProf) {
        EXPECT_EQUAL(fixture.unReducedData[test][fixture.starts[iProf]],
                     fixture.data[indexer.at(iProf)]);
      }
    }
    //// reducing copy constructor
    for (size_t iSelect = 0; iSelect < fixture.selections.size();
        ++iSelect) {
      SECTION("reduced indexing of profile " + indexingType
          + " indices with " + fixture.selectionTypes[iSelect] + " selected") {
        std::vector<bool> select = fixture.selections[iSelect];
        std::vector<size_t> reducedData =
          fixture.ReducedData[test][iSelect];

        DataIndexer reducedIndexer(indexer, select);
        std::cout << "reducedIndexer indices: " << std::endl;
        for (size_t iProf = 0; iProf < reducedIndexer.n_obs(); ++iProf) {
          std::cout << iProf << ": ";
          for (size_t iLevel = 0; iLevel < reducedIndexer.length(iProf);
              ++iLevel) {
            std::cout << reducedIndexer.at(iProf, iLevel) << " ";
          }
          std::cout << std::endl;
        }
        EXPECT_EQUAL(reducedIndexer.n_obs(),
                     fixture.nProfiles[test][iSelect]);
        EXPECT_EQUAL(reducedIndexer.n_levels(),
                     fixture.nLevels[test][iSelect]);
        {
          size_t iProfReduced = 0;
          size_t iLocReduced = 0;
          for (size_t iProf = 0; iProf < indexer.n_obs(); ++iProf) {
            size_t iLevelReduced = 0;
            for (size_t iLevel = 0; iLevel < indexer.length(iProf); ++iLevel) {
              if (select[indexer(iProf, iLevel)]) {
                EXPECT_EQUAL(reducedData[iLocReduced],
                    fixture.data[reducedIndexer(iProfReduced, iLevelReduced)]);
                iLevelReduced++;
                iLocReduced++;
              }
            }
            if (iLevelReduced != 0) iProfReduced++;
          }
        }
      }
    }
  }
}

CASE("Test Nemo Feedback surface data") {
  // Constrcted using deep copy of the data
  std::vector<size_t> indices = fixture.indices.at(0);
  auto indexer = std::make_shared<DataIndexer>(indices,
      fixture.unReducedData[0].size());
  Data<size_t>      feedbackData(indexer, fixture.data);
  Data<float>       feedbackFloatData(indexer, fixture.floatData);
  Data<double>      feedbackDoubleData(indexer, fixture.doubleData);
  Data<int32_t>     feedbackIntData(indexer, fixture.intData);
  Data<std::string> feedbackStringData(indexer, fixture.stringData);

  SECTION("Feedback surface data editing") {
    for (size_t iLoc = 0; iLoc < feedbackData.n_obs(); ++iLoc) {
      auto initialValue = fixture.data[iLoc];
      feedbackData[iLoc] = 42;
      auto result = feedbackData[iLoc];
      EXPECT_EQUAL(result, 42);
      EXPECT_EQUAL(initialValue, fixture.data[iLoc]);
    }
  }

  SECTION("Feedback surface string data editing") {
    for (size_t iLoc = 0; iLoc < feedbackStringData.n_obs(); ++iLoc) {
      std::string initialValue = fixture.stringData[iLoc];
      feedbackStringData[iLoc] = "X";
      std::string resultString = feedbackStringData[iLoc];
      EXPECT_EQUAL(resultString, "X");
      resultString = fixture.stringData[iLoc];
      EXPECT_EQUAL(resultString, initialValue);
    }
  }

  SECTION("Feedback surface raw data retrieval") {
    auto values = feedbackData.raw();
    for (size_t iLoc = 0; iLoc < values.size(); ++iLoc) {
      EXPECT_EQUAL(values.at(iLoc), fixture.data.at(iLoc));
    }
  }

  SECTION("Feedback surface shallow constructor") {
    auto dataPtr = std::make_shared<std::vector<size_t>>(fixture.data);
    Data<size_t> feedbackDataShallow(indexer, dataPtr);
    // editing shallow copy
    for (size_t iLoc = 0; iLoc < feedbackDataShallow.n_obs(); ++iLoc) {
      feedbackDataShallow[iLoc] = 42;
      EXPECT_EQUAL((*dataPtr)[iLoc], 42);
    }
  }

  SECTION("Feedback surface editing deep copy doesn't affect other copies") {
    auto dataPtr = std::make_shared<std::vector<size_t>>(fixture.data);
    Data<size_t> feedbackDataShallow(indexer, dataPtr);
    Data<size_t> feedbackDataDeep = feedbackDataShallow.deep_copy();
    for (size_t iLoc = 0; iLoc < feedbackDataShallow.n_obs(); ++iLoc) {
      feedbackDataShallow[iLoc] = 42;
      EXPECT_EQUAL(feedbackDataDeep[iLoc], fixture.data[iLoc]);
    }
  }

  for (size_t test = 0; test < fixture.indices.size(); ++test) {
    std::string indexingType = fixture.indexingTypes.at(test);
    std::vector<size_t> indices = fixture.indices.at(test);
    // Deep copy of the data
    auto indexer = std::make_shared<DataIndexer>(indices,
        fixture.unReducedData[0].size());
    Data<size_t>      feedbackData(indexer, fixture.data);
    Data<float>       feedbackFloatData(indexer, fixture.floatData);
    Data<double>      feedbackDoubleData(indexer, fixture.doubleData);
    Data<int32_t>     feedbackIntData(indexer, fixture.intData);
    Data<std::string> feedbackStringData(indexer, fixture.stringData);
    for (size_t iSelect = 0; iSelect < fixture.selections.size();
        ++iSelect) {
      SECTION("update surface indexing " + indexingType
          + " indices with " + fixture.selectionTypes[iSelect] + " selected") {
        std::vector<bool> select = fixture.selections[iSelect];
        // update the float data using the reducing constructor for this test
        Data<float> feedbackReducedFloatData(feedbackFloatData, select);
        // update the indexer to point at the floatData indexer
        feedbackData.update_indexer(feedbackReducedFloatData);
        std::vector<size_t> reducedData =
          fixture.ReducedData[test][iSelect];
        EXPECT_EQUAL(feedbackData.n_obs(), reducedData.size());
        size_t nLevelsReduced = reducedData.size() > 0 ? 1 : 0;
        EXPECT_EQUAL(feedbackData.n_levels(), nLevelsReduced);
        for (size_t iLoc = 0; iLoc < feedbackData.n_obs(); ++iLoc) {
            EXPECT_EQUAL(feedbackData.at(iLoc), reducedData.at(iLoc));
        }
        // update the indexer to point to the original indexer
        feedbackData.update_indexer(indexer);
      }
    }
  }
}

CASE("Test Nemo Feedback profile data") {
  // Constrcted using deep copy of the data
  std::vector<size_t> indices = fixture.indices.at(0);
  auto indexer = std::make_shared<DataIndexer>(fixture.unreducedNProfiles,
      fixture.unreducedNLevels, fixture.unReducedData[0].size(), fixture.starts,
      indices);
  Data<size_t>      feedbackData(indexer, fixture.data);
  Data<float>       feedbackFloatData(indexer, fixture.floatData);
  Data<double>      feedbackDoubleData(indexer, fixture.doubleData);
  Data<int32_t>     feedbackIntData(indexer, fixture.intData);
  Data<std::string> feedbackStringData(indexer, fixture.stringData);

  SECTION("Feedback profile data editing") {
      size_t iLoc = 0;
      for (size_t iProf = 0; iProf < feedbackData.n_obs(); ++iProf) {
        for (size_t iLevel = 0; iLevel < feedbackData.length(iProf); ++iLevel) {
          auto initialValue = fixture.data.at(iLoc);
          feedbackData(iProf, iLevel) = 42;
          auto result = feedbackData.at(iProf, iLevel);
          EXPECT_EQUAL(result, 42);
          EXPECT_EQUAL(initialValue, fixture.data.at(iLoc));
          iLoc++;
        }
      }
    }

  SECTION("Feedback profile string data editing") {
    size_t iLoc = 0;
    for (size_t iProf = 0; iProf < feedbackStringData.n_obs(); ++iProf) {
      for (size_t iLevel = 0; iLevel < feedbackStringData.length(iProf);
          ++iLevel) {
        auto initialValue = fixture.stringData[iLoc];
        feedbackStringData(iProf, iLevel) = "X";
        auto result = feedbackStringData(iProf, iLevel);
        EXPECT_EQUAL(result, "X");
        EXPECT_EQUAL(initialValue, fixture.stringData[iLoc]);
        iLoc++;
      }
    }
  }

  SECTION("Feedback profile raw data retrieval") {
    auto values = feedbackData.raw();
    for (size_t iLoc = 0; iLoc < values.size(); ++iLoc) {
      EXPECT_EQUAL(values.at(iLoc), fixture.data.at(iLoc));
    }
  }

  SECTION("Feedback profile raw profile data retrieval") {
    size_t iLoc = 0;
    for (size_t iProf = 0; iProf < feedbackData.n_obs(); ++iProf) {
      auto values = feedbackData.raw_profile(iProf);
      for (size_t iLevel = 0; iLevel < feedbackData.length(iProf); ++iLevel) {
        EXPECT_EQUAL(values.at(iLevel), fixture.data.at(iLoc));
        iLoc++;
      }
    }
  }

  SECTION("Feedback profile shallow constructor") {
    auto dataPtr = std::make_shared<std::vector<size_t>>(fixture.data);
    Data<size_t> feedbackDataShallow(indexer, dataPtr);
    size_t iLoc = 0;
    for (size_t iProf = 0; iProf < feedbackDataShallow.n_obs(); ++iProf) {
      for (size_t iLevel = 0; iLevel < feedbackDataShallow.length(iProf);
          ++iLevel) {
        feedbackDataShallow(iProf, iLevel) = 42;
        EXPECT_EQUAL((*dataPtr)[iLoc], 42);
        iLoc++;
      }
    }
  }

  SECTION("Feedback profile editing deep copy doesn't affect other copies") {
    auto dataPtr = std::make_shared<std::vector<size_t>>(fixture.data);
    Data<size_t> feedbackDataShallow(indexer, dataPtr);
    Data<size_t> feedbackDataDeep = feedbackDataShallow.deep_copy();
    EXPECT_EQUAL(feedbackDataShallow.n_obs(), feedbackDataDeep.n_obs());
    EXPECT_EQUAL(feedbackDataShallow.n_levels(), feedbackDataDeep.n_levels());
    size_t iLoc = 0;
    for (size_t iProf = 0; iProf < feedbackDataShallow.n_obs(); ++iProf) {
      for (size_t iLevel = 0; iLevel < feedbackDataShallow.length(iProf);
          ++iLevel) {
        feedbackDataShallow(iProf, iLevel) = 42;
        EXPECT_EQUAL(feedbackDataDeep.at(iProf, iLevel), fixture.data[iLoc]);
        iLoc++;
      }
    }
  }

  for (size_t test = 0; test < fixture.indices.size(); ++test) {
    std::string indexingType = fixture.indexingTypes.at(test);
    std::vector<size_t> indices = fixture.indices.at(test);
    // Deep copy of the data
    auto indexer = std::make_shared<DataIndexer>(fixture.unreducedNProfiles,
        fixture.unreducedNLevels, fixture.unReducedData[0].size(),
        fixture.starts, indices);
    Data<size_t>      feedbackData(indexer, fixture.data);
    Data<float>       feedbackFloatData(indexer, fixture.floatData);
    Data<double>      feedbackDoubleData(indexer, fixture.doubleData);
    Data<int32_t>     feedbackIntData(indexer, fixture.intData);
    Data<std::string> feedbackStringData(indexer, fixture.stringData);
    for (size_t iSelect = 0; iSelect < fixture.selections.size();
        ++iSelect) {
      SECTION("update profile indexing " + indexingType + " indices with "
          + fixture.selectionTypes[iSelect] + " selected") {
        std::vector<bool> select = fixture.selections[iSelect];
        // update the float data using the reducing constructor
        Data<float> feedbackReducedFloatData(feedbackFloatData, select);
        // update the indexer to point at the floatData indexer
        feedbackData.update_indexer(feedbackReducedFloatData);
        std::vector<size_t> reducedData =
          fixture.ReducedData[test][iSelect];
        EXPECT_EQUAL(
            feedbackData.n_obs(), fixture.nProfiles[test][iSelect]);
        EXPECT_EQUAL(
            feedbackData.n_levels(), fixture.nLevels[test][iSelect]);
        size_t iLoc = 0;
        for (size_t iProf = 0; iProf < feedbackData.n_obs(); ++iProf) {
          for (size_t iLevel = 0; iLevel < feedbackData.length(iProf);
              ++iLevel) {
            EXPECT_EQUAL(feedbackData.at(iProf, iLevel), reducedData.at(iLoc));
            iLoc++;
          }
        }
        // update the indexer to point to the original indexer
        feedbackData.update_indexer(indexer);
      }
    }
  }
}

}  // namespace test
}  // namespace feedback_io
}  // namespace nemo_feedback

int main(int argc, char** argv) {
    return eckit::testing::run_tests(argc, argv);
}
