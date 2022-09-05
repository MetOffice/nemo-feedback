/*
 * (C) British Crown Copyright 2022 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include <vector>
#include <algorithm>
#include <numeric>
#include <type_traits>

#include "eckit/exception/Exceptions.h"
#include "oops/util/missingValues.h"

#include "nemo-feedback/NemoFeedbackReduce.h"
#include "nemo-feedback/NemoFeedbackWriter.h"

namespace nemo_feedback {

NemoFeedbackReduce::NemoFeedbackReduce(const size_t n_obs, const size_t n_obs_to_write,
    const std::vector<bool>& to_write,
    const std::vector<size_t> & record_starts,
    const std::vector<size_t> & record_counts) :
    n_obs_(n_obs), n_obs_to_write_(n_obs_to_write), to_write_(to_write),
    unreduced_starts(record_starts), unreduced_counts(record_counts) {
  reduced_starts.reserve(n_obs_);
  reduced_counts.reserve(n_obs_);
  for (int i = 0; i < n_obs_; ++i) {
    size_t reclen = 0;
    for (int l = record_starts[i]; l < record_starts[i]+record_counts[i]; ++l) {
      if (to_write_[l]) {
        if (i >= reduced_starts.size()) {
          if (i == 0) {
            reduced_starts.push_back(l);
          } else {
            reduced_starts.push_back(
                reduced_starts[i-1]+reduced_counts[i-1]);
          }
        }
        reclen++;
      }
    }
    reduced_counts.push_back(reclen);
  }
}

template <typename T>
std::vector<T> NemoFeedbackReduce::reduce_data(
    const std::vector<T> & data_in,
    const bool change_fillvalues) {
  // profile data, n_obs dimension var - using the n_obs should be fine
  // surface data, n_locs/n_obs dimension var, n_locs = n_obs
  // so using n_obs should be fine
  auto missing_value = util::missingValue(T(0));
  std::vector<T> data_out(n_obs_to_write_);
  int j = 0;
  for (int i = 0; i < n_obs_; ++i) {
    if (to_write_[i]) {
      data_out[j++] = data_in[i];
      if (change_fillvalues && (data_in[i] == missing_value)){
        data_out[j-1] = static_cast<T>(NemoFeedbackWriter::double_fillvalue);
      }
    }
  }
  return data_out;
}

template std::vector<int> NemoFeedbackReduce::reduce_data<int>(
    const std::vector<int> & data_in, const bool change_fillvalues);
template std::vector<float> NemoFeedbackReduce::reduce_data<float>(
    const std::vector<float> & data_in, const bool change_fillvalues);
template std::vector<double> NemoFeedbackReduce::reduce_data<double>(
    const std::vector<double> & data_in, const bool change_fillvalues);

std::vector<std::string> NemoFeedbackReduce::reduce_data(
    const std::vector<std::string> & data_in) {
  std::vector<std::string> data_out(n_obs_to_write_);
  int j = 0;
  for (int i = 0; i < n_obs_; ++i) {
    if (to_write_[i]) {
      data_out[j++] = data_in[i];
    }
  }
  return data_out;
}

template <typename T>
void NemoFeedbackReduce::reduce_profile_data(
    const std::vector<T> & data_in,
    std::vector<T> & data_out,
    const bool change_fillvalues
    ) {
  // check that the size of the input data matches the unreduced_counts
  const size_t n_unred_prof_obs = std::accumulate(unreduced_counts.begin(), unreduced_counts.end(),
      decltype(reduced_counts)::value_type(0));
  if (data_in.size() != n_unred_prof_obs) {
        throw eckit::BadValue("NemoFeedbackReduce:: bad counts or input data size: "
            + std::to_string(data_in.size()) + " with nLocs " + std::to_string(n_unred_prof_obs),
            Here());
  }
  // with profile data n_obs != n_locs, and so we setup new record_starts and
  // counts based on the new data vector.
  data_out.clear();
  const size_t n_prof_obs = std::accumulate(reduced_counts.begin(), reduced_counts.end(),
      decltype(reduced_counts)::value_type(0));
  data_out.reserve(n_prof_obs);
  auto missing_value = util::missingValue(T(0));
  for (int iprof = 0; iprof < n_obs_; ++iprof) {
    size_t reclen = 0;
    for (int l = unreduced_starts[iprof]; l < unreduced_starts[iprof]+unreduced_counts[iprof]; ++l) {
      if (to_write_[l]) {
        if (reclen++ == reduced_counts[iprof]) break;
        data_out.push_back(data_in[l]);
        if (change_fillvalues && (data_in[l] == missing_value))
             data_out.back() = static_cast<T>(NemoFeedbackWriter::double_fillvalue);
      }
    }
  }
}

template void NemoFeedbackReduce::reduce_profile_data(
    const std::vector<int> & data_in,
    std::vector<int> & data_out,
    const bool change_fillvalues
    );
template void NemoFeedbackReduce::reduce_profile_data(
    const std::vector<float> & data_in,
    std::vector<float> & data_out,
    const bool change_fillvalues
    );
template void NemoFeedbackReduce::reduce_profile_data(
    const std::vector<double> & data_in,
    std::vector<double> & data_out,
    const bool change_fillvalues
    );

}   // namespace nemo_feedback
