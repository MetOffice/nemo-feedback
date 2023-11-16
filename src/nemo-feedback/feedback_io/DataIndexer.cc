/*
 * (C) British Crown Copyright 2023 Met Office
 */


#include "nemo-feedback/feedback_io/DataIndexer.h"

#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <sstream>

#include "eckit/exception/Exceptions.h"

namespace nemo_feedback {
namespace feedback_io {

DataIndexer::DataIndexer() : nObs_(0), nLevels_(0), sourceDataSize_(0),
    indices_(std::vector<size_t>()), starts_(std::vector<size_t>()),
    counts_(std::vector<size_t>()) {
      validate();
}

DataIndexer::DataIndexer(const size_t nObs, const size_t nLevels,
                         const size_t sourceDataSize,
                         const std::vector<size_t> starts,
                         const std::vector<size_t> indices) :
  nObs_(nObs), nLevels_(nLevels),  sourceDataSize_(sourceDataSize),
  indices_(indices), starts_(starts) {
  counts_.reserve(nObs);
  for (size_t iProf = 0; iProf < nObs - 1; ++iProf) {
    counts_.push_back(starts[iProf +1] - starts[iProf]);
  }
  // cannot infer the number in the final profile from the start vector, must
  // use the indices vector.
  counts_.push_back(indices.size() - starts[nObs - 1]);
  validate();
}

DataIndexer::DataIndexer(const std::vector<size_t> indices,
                         const size_t sourceDataSize) :
    nObs_(indices.size()), nLevels_(1), sourceDataSize_(sourceDataSize),
    indices_(indices), starts_(indices.size()), counts_(indices.size(), 1) {
  std::iota(starts_.begin(), starts_.end(), 0);
  validate();
}

DataIndexer::DataIndexer(const DataIndexer& superset,
  const std::vector<bool> select) : nObs_(0), nLevels_(0),
    sourceDataSize_(superset.sourceDataSize_),
    indices_(), starts_(), counts_() {
  ASSERT_MSG(superset.indices_.size() == select.size(),
      DataIndexer::className()
      + ": source indices and selection vector size mismatch.");
  counts_.reserve(superset.nObs_);
  starts_.reserve(superset.nObs_);
  size_t iProfReduced = 0;
  for (size_t iProf = 0; iProf < superset.nObs_; ++iProf) {
    size_t nValidLevels = 0;
    for (size_t iLevel = 0; iLevel < superset.length(iProf); ++iLevel) {
      if (select[superset(iProf, iLevel)]) {
        nValidLevels++;
      }
    }
    if (nValidLevels > 0) {
      // setup nObs and nLevels
      nObs_++;
      if (nValidLevels > nLevels_) nLevels_ = nValidLevels;
      // setup counts
      counts_.push_back(nValidLevels);
      // setup indices
      for (size_t iLevel = 0; iLevel < superset.length(iProf); ++iLevel) {
        size_t srcIndex = superset(iProf, iLevel);
        if (select[srcIndex]) {
          indices_.push_back(srcIndex);
        }
      }
      // setup starts
      if (iProfReduced == 0) {
        starts_.push_back(0);
      } else {
        starts_.push_back(starts_[iProfReduced-1] + counts_[iProfReduced-1]);
      }
      iProfReduced++;
    }
  }

  // sanity check that the amount of data has reduced
  ASSERT_MSG(superset.nObs_ >= nObs_,
      DataIndexer::className() + ": superset smaller than subset "
      + std::to_string(superset.nObs_) + " !>= " + std::to_string(nObs_));
  ASSERT_MSG(superset.nLevels_ >= nLevels_,
      DataIndexer::className() + ": superset has fewer levels than subset "
      + std::to_string(superset.nLevels_) + " !>= " + std::to_string(nLevels_));

  // sanity check that our new indexer is a valid indexer
  validate();
}

void DataIndexer::validate() {
  std::ostringstream message;
  message << DataIndexer::className() << ": starts vector size mismatch: "
          << starts_.size() << " == " << nObs_;
  ASSERT_MSG(starts_.size() == nObs_, message.str());
  message.str(std::string());

  message << DataIndexer::className() << ": counts vector size mismatch: "
          << counts_.size() << " == " << nObs_;
  ASSERT_MSG(counts_.size() == nObs_, message.str());
  message.str(std::string());

  size_t maxLevel = 0;
  if (nObs_ != 0) {
    maxLevel = (*max_element(counts_.begin(), counts_.end()));
  }
  message << DataIndexer::className()
          << ": max profile size not equal to nLevels: "
          << maxLevel << " == " << nLevels_;
  ASSERT_MSG(maxLevel == nLevels_, message.str());
  message.str(std::string());

  message << DataIndexer::className()
          << ": indices vector mismatch with nObs, nLevels: " << indices_.size()
          << " > " << nObs_ << " * " << nLevels_;
  ASSERT_MSG(indices_.size() <= nObs_ * nLevels_, message.str());
  message.str(std::string());

  const size_t nLocations = std::accumulate(counts_.begin(), counts_.end(),
      size_t(0));
  message << DataIndexer::className()
          << ": indices vector size mismatch with counts: " << indices_.size()
          << " != " << nLocations;
  ASSERT_MSG(indices_.size() == nLocations, message.str());
  message.str(std::string());

  if ((nObs_ > 0) && (nLevels_ > 0)) {
    message << DataIndexer::className()
            << ": last element does not fit in indices vector: "
            << indices_.size() << " < " << starts_[nObs_-1] + counts_[nObs_-1];
    ASSERT_MSG(indices_.size() == starts_[nObs_-1] + counts_[nObs_-1],
        message.str());
    message.str(std::string());
  }
}
}  // namespace feedback_io
}  // namespace nemo_feedback
