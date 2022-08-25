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
      const std::vector<bool>& to_write) :
    n_obs_(n_obs), n_obs_to_write_(n_obs_to_write), to_write_(to_write) {};

  /// \brief remove unwanted data according to the `to_write_` vector
  template <typename T>
  std::vector<T> reduce_data(
      const std::vector<T> & data_in);

  /// \brief remove unwanted data according to the `to_write_` vector for
  ///        profile data
  template <typename T>
  void reduce_profile_data(
    const std::vector<size_t> & record_starts,
    const std::vector<size_t> & record_counts,
    const std::vector<T> & data_in,
    std::vector<size_t> & record_starts_out,
    std::vector<size_t> & record_counts_out,
    std::vector<T> & data_out);

  const size_t n_obs_;
  const size_t n_obs_to_write_;
  const std::vector<bool> to_write_;
};

}
