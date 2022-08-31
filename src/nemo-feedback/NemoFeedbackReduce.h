/*
 * (C) British Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#pragma once

#include <vector>

namespace nemo_feedback {

class NemoFeedbackReduce {
 public:
  NemoFeedbackReduce(const size_t n_obs, const size_t n_obs_to_write,
      const std::vector<bool>& to_write,
      const std::vector<size_t> & record_starts,
      const std::vector<size_t> & record_counts);

  /// \brief remove unwanted data according to the `to_write_` vector
  template <typename T>
  std::vector<T> reduce_data(
      const std::vector<T> & data_in);

  /// \brief remove unwanted data according to the `to_write_` vector for
  ///        profile data
  template <typename T>
  void reduce_profile_data(
      const std::vector<T> & data_in,
      std::vector<T> & data_out,
      const bool change_fillvalues = true);

  /// \brief remove unwanted data according to the `validObs` vector
  template <typename T>
  std::vector<T> reduce_via_accessor(
      const std::vector<T> & data_in,
      const std::vector<size_t> & validObs,
    const bool change_fillvalues = true);

  const size_t n_obs_;
  const size_t n_obs_to_write_;
  const std::vector<bool> to_write_;
  std::vector<size_t> reduced_starts;
  std::vector<size_t> reduced_counts;
  std::vector<size_t> unreduced_starts;
  std::vector<size_t> unreduced_counts;
};

}
