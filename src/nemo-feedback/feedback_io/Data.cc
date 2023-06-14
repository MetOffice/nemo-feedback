/*
 * (C) British Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include <iostream>
#include <sstream>
#include <string>
#include <memory>
#include <numeric>
#include <utility>

#include "eckit/exception/Exceptions.h"

#include "nemo-feedback/feedback_io/Data.h"

namespace nemo_feedback {
namespace feedback_io {

template <class T>
void Data<T>::validate() const {
  // Note, this check assumes each datum is indexed only once
  std::ostringstream message;
  message << Data::className() << ": indexer too small for data: "
          << indexer_->n_locations() << " > " << data_->size();
  ASSERT_MSG(indexer_->n_locations() <= data_->size(), message.str());
  message.str(std::string());
  message << Data::className() << ": indexer created for a different vector: "
          << indexer_->n_source_data() << " != " << data_->size();
  ASSERT_MSG(indexer_->n_source_data() == data_->size(), message.str());
  message.str(std::string());
}

template <class T>
Data<T> Data<T>::deep_copy() const {
  std::vector<T> data = this->raw();
  std::vector<size_t> indices(data.size());
  std::iota(indices.begin(), indices.end(), size_t(0));

  std::vector<size_t> starts;
  starts.push_back(0);
  for (size_t iProfile = 1; iProfile < this->indexer_->n_obs(); ++iProfile) {
    starts.push_back(starts[iProfile-1] + this->indexer_->length(iProfile-1));
  }

  DataIndexer indexer(
      this->indexer_->n_obs(), this->indexer_->n_levels(), data.size(),
      starts, indices);

  return Data(std::move(indexer), std::move(data));
}

template <class T>
std::vector<T> Data<T>::raw() const {
  std::vector<T> result;
  result.reserve(data_->size());
  for (size_t iProfile = 0; iProfile < this->n_obs(); ++iProfile) {
    for (size_t iLevel = 0; iLevel < this->length(iProfile); ++iLevel) {
      result.emplace_back((*this)(iProfile, iLevel));
    }
  }
  return result;
}

template <class T>
std::vector<T> Data<T>::raw_surface() const {
  std::vector<T> result;
  result.reserve(this->indexer_->n_obs());
  for (size_t iLoc = 0; iLoc < this->indexer_->n_obs(); ++iLoc) {
    result.emplace_back((*this)[iLoc]);
  }
  return result;
}

template <class T>
std::vector<T> Data<T>::raw_profile(size_t iProfile) const {
  std::vector<T> result;
  result.reserve(this->length(iProfile));
  for (size_t iLevel = 0; iLevel < this->length(iProfile); ++iLevel) {
    result.emplace_back((*this)(iProfile, iLevel));
  }
  return result;
}

template class Data<std::string>;
template class Data<size_t>;
template class Data<int32_t>;
template class Data<float>;
template class Data<double>;
}  // namespace feedback_io
}  // namespace nemo_feedback
