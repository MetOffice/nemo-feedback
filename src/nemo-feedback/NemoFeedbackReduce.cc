/*
 * (C) British Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include <vector>
#include <algorithm>

#include "nemo-feedback/NemoFeedbackReduce.h"

namespace nemo_feedback {

template <typename T>
std::vector<T> NemoFeedbackReduce::reduce_data(
    const std::vector<T> & data_in) {
  // profile data, n_obs dimension var - using the n_obs should be fine
  // surface data, n_locs/n_obs dimension var, n_locs = n_obs
  // so using n_obs should be fine
  std::vector<T> data_out(n_obs_to_write_);
  int j = 0;
  for (int i = 0; i < n_obs_; ++i) {
    if (to_write_[i]) {
      data_out[j++] = data_in[i];
    }
  }
  return data_out;
}

template std::vector<int> NemoFeedbackReduce::reduce_data<int>(
    const std::vector<int> & data_in);
template std::vector<float> NemoFeedbackReduce::reduce_data<float>(
    const std::vector<float> & data_in);
template std::vector<double> NemoFeedbackReduce::reduce_data<double>(
    const std::vector<double> & data_in);

template <typename T>
void NemoFeedbackReduce::reduce_profile_data(
    const std::vector<size_t> & record_starts,
    const std::vector<size_t> & record_counts,
    const std::vector<T> & data_in,
    std::vector<size_t> & record_starts_out,
    std::vector<size_t> & record_counts_out,
    std::vector<T> & data_out
    ) {
  // with profile data n_obs != n_locs, and so we setup new record_starts and
  // counts based on the new data vector.
  data_out.reserve((std::count(to_write_.begin(), to_write_.end(), true)));
  record_starts_out.reserve(n_obs_);
  record_counts_out.reserve(n_obs_);
  for (int i = 0; i < n_obs_; ++i) {
    size_t reclen = 0;
    for (int l = record_starts[i]; l < record_starts[i]+record_counts[i]; ++l) {
      if (to_write_[l]) {
        if (i >= record_starts_out.size()) {
          if (i == 0) {
            record_starts_out.push_back(l);
          } else {
            record_starts_out.push_back(
                record_starts_out[i-1]+record_counts_out[i-1]);
          }
        }
        data_out.push_back(data_in[l]);
        reclen++;
      }
    }
    record_counts_out.push_back(reclen);
  }
}

}
