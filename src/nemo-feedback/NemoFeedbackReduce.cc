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
#include "oops/util/Logger.h"

#include "nemo-feedback/NemoFeedbackReduce.h"
#include "nemo-feedback/NemoFeedbackWriter.h"

namespace nemo_feedback {

NemoFeedbackReduce::NemoFeedbackReduce(const size_t n_obs,
    const size_t n_surf_obs_to_write,
    const std::vector<bool>& to_write,
    const std::vector<size_t> & record_starts,
    const std::vector<size_t> & record_counts) :
    n_obs_(n_obs), n_surf_obs_to_write_(n_surf_obs_to_write), to_write_(to_write),
    unreduced_starts(record_starts), unreduced_counts(record_counts) {
  oops::Log::trace() << "NemoFeedbackReduce constructor" << std::endl;
  reduced_starts.reserve(n_obs_);
  reduced_counts.reserve(n_obs_);
  for (int i = 0; i < n_obs_; ++i) {
    size_t reclen = 0;
    for (int l = record_starts[i]; l < record_starts[i]+record_counts[i]; ++l) {
      if (to_write_[l])
        reclen++;
    }
    reduced_counts.push_back(reclen);
  }
  reduced_starts.push_back(0);
  for (int i = 1; i < n_obs; ++i) {
    reduced_starts.push_back(reduced_starts[i-1] + reduced_counts[i-1]);
  }
  const size_t total_reduced_counts = std::accumulate(reduced_counts.begin(),
      reduced_counts.end(), decltype(reduced_counts)::value_type(0));
  const size_t total_to_write = std::count(to_write.begin(), to_write.end(), true);
  if ( total_reduced_counts != total_to_write ) {
      std::ostringstream err_stream;
      err_stream << "nemo_feedback::NemoFeedbackReduce::constructor "
                 << "total reduced counts doesn't match total number of obs to write"
                 << total_reduced_counts << " != " << total_to_write;
      throw eckit::BadValue(err_stream.str(), Here());
  }
}

template <typename T>
std::vector<T> NemoFeedbackReduce::reduce_data(
    const std::vector<T> & data_in,
    const bool change_fillvalues) {
  oops::Log::trace() << "NemoFeedbackReduce::reduce_data" << std::endl;
  // profile data, n_obs dimension var - using the n_obs should be fine
  // surface data, n_locs/n_obs dimension var, n_locs = n_obs
  // so using n_obs should be fine
  auto missing_value = util::missingValue(T(0));
  std::vector<T> data_out(n_surf_obs_to_write_);
  int j = 0;
  for (int i = 0; i < n_obs_; ++i) {
    if (to_write_[i]) {
      data_out[j++] = data_in[i];
      if (change_fillvalues && (data_in[i] == missing_value)) {
        data_out[j-1] = static_cast<T>(typeToFill::value<T>());
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
  std::vector<std::string> data_out(n_surf_obs_to_write_);
  oops::Log::trace() << "NemoFeedbackReduce::reduce_data<string>" << std::endl;
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
  oops::Log::trace() << "NemoFeedbackReduce::reduce_profile_data" << std::endl;
  // check that the size of the input data is large enough compared to the
  // amount of data expected in unreduced_counts. The data described in
  // unreduced_counts may be smaller than the size of the input data because
  // unreduced_counts does not include data from whole profiles that are
  // not being written to file.
  const size_t n_unred_prof_locs = std::accumulate(unreduced_counts.begin(),
      unreduced_counts.end(), decltype(reduced_counts)::value_type(0));
  if (data_in.size() < n_unred_prof_locs) {
        throw eckit::BadValue(
            "NemoFeedbackReduce:: bad counts or input data size: "
            + std::to_string(data_in.size()) + " with nLocs "
            + std::to_string(n_unred_prof_locs), Here());
  }
  // with profile data n_obs is the number of profiles and hence is not equal
  // to n_locs, which is the number of data points, and so we set up new
  // record_starts and counts based on the new data vector.
  data_out.clear();
  const size_t n_prof_locs = std::accumulate(reduced_counts.begin(),
      reduced_counts.end(), decltype(reduced_counts)::value_type(0));
  data_out.reserve(n_prof_locs);
  auto missing_value = util::missingValue(T(0));
  for (int iprof = 0; iprof < n_obs_; ++iprof) {
    size_t reclen = 0;
    size_t unreduced_recend = unreduced_starts[iprof]+unreduced_counts[iprof];
    for (size_t l = unreduced_starts[iprof]; l < unreduced_recend; ++l) {
      if (to_write_[l]) {
        if (reclen++ == reduced_counts[iprof]) break;
        data_out.push_back(data_in[l]);
        if (change_fillvalues && (data_in[l] == missing_value)) {
          data_out.back() = static_cast<T>(typeToFill::value<T>());
        }
      }
    }
  }
  if (reduced_counts.at(n_obs_-1) + reduced_starts.at(n_obs_-1) != data_out.size()) {
      std::ostringstream err_stream;
      err_stream << "nemo_feedback::NemoFeedbackReduce::reduce_profile_data "
                 << "index range out of bounds n_prof_locs "
                 << n_obs_ << " reduced_starts.size() "
                 << reduced_starts.size() << " counts: "
                 << reduced_counts.at(n_obs_-1) << " + "
                 << reduced_starts.at(n_obs_-1) << " = "
                 << reduced_counts.at(n_obs_-1) +
                    reduced_starts.at(n_obs_-1)
                 << " >= " << data_out.size();
      throw eckit::BadValue(err_stream.str(), Here());
  }
}

template void NemoFeedbackReduce::reduce_profile_data(
    const std::vector<int> & data_in,
    std::vector<int> & data_out,
    const bool change_fillvalues);
template void NemoFeedbackReduce::reduce_profile_data(
    const std::vector<float> & data_in,
    std::vector<float> & data_out,
    const bool change_fillvalues);
template void NemoFeedbackReduce::reduce_profile_data(
    const std::vector<double> & data_in,
    std::vector<double> & data_out,
    const bool change_fillvalues);

}   // namespace nemo_feedback
