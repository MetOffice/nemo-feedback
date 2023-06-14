/*
 * (C) British Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#pragma once

#include <vector>
#include <string>

namespace nemo_feedback {
namespace feedback_io {

  /// \brief Indexing of a 1D vector of observation data by the sorted
  ///        (profile, level) or (location) for surface information
class DataIndexer {
 public:
  static std::string className()
    { return "nemo_feedback::feedback_io::DataIndexer"; }
  /// \brief empty constructor
  DataIndexer();
  /// \brief profile data indexing constructor
  DataIndexer(const size_t nObs, const size_t nLevels,
              const size_t sourceDataSize,
              const std::vector<size_t> starts,
              const std::vector<size_t> indices);
  /// \brief surface data indexing constructor
  explicit DataIndexer(const std::vector<size_t> indices,
                       const size_t sourceDataSize);
  /// \brief force a default copy constructor
  DataIndexer(const DataIndexer&) = default;
  /// \brief Copy constructor for a new indexer which is a subset of another
  /// indexer
  DataIndexer(const DataIndexer& superset, const std::vector<bool> select);
  size_t operator()(size_t iLocation) const
    { return indices_[starts_[iLocation]]; }
  size_t operator()(size_t iProfile, size_t iLevel) const
    { return indices_[starts_[iProfile] + iLevel]; }
  size_t at(size_t iLocation) const
    { return indices_.at(starts_.at(iLocation)); }
  size_t at(size_t iProfile, size_t iLevel) const
    { return indices_.at(starts_.at(iProfile) + iLevel); }
  size_t length(size_t iProfile) const { return counts_.at(iProfile); }
  size_t n_obs() const { return nObs_; }
  size_t n_levels() const { return nLevels_; }
  size_t n_locations() const { return indices_.size(); }
  size_t n_source_data() const { return sourceDataSize_; }

 private:
  void validate();
  size_t nObs_;
  size_t nLevels_;
  const size_t sourceDataSize_;
  std::vector<size_t> indices_;
  std::vector<size_t> starts_;
  std::vector<size_t> counts_;
};
}  // namespace feedback_io
}  // namespace nemo_feedback
