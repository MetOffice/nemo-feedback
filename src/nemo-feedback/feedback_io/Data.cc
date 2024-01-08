/*
 * (C) British Crown Copyright 2024 Met Office
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
          << indexer_->n_source_data() << " < " << data_->size();
  ASSERT_MSG(indexer_->n_source_data() >= data_->size(), message.str());
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
template class Data<feedback_io::QC::Level>;
template class Data<float>;
template class Data<double>;

/// \brief update the whole report QC information for a profile based on the
///        current whole variable QC information
feedback_io::QC::Level wholeReportUpdate(feedback_io::QC::Level current, const
    size_t length, const size_t nGoodObs, const size_t nBadObs, const size_t
    nBadDoNotAssimilate, const size_t nDoNotAssimilate) {
  // If good obs were previously found, the report is good
  if (current == QC::Level::Good) {
    return QC::Level::Good;
  }
  // Upgrade any reports if they contain any good obs (goodDoNotAssimilate ->
  // good)
  if (nGoodObs > 0) {
    return QC::Level::Good;
  }
  // Check currently bad, do-not-assimilate or not-yet-checked reports and
  // downgrade where necessary
  if (nBadObs == length) {
    return QC::Level::Bad;
  } else if (nBadDoNotAssimilate == length) {
    // No good observations, all bad-do-not-assimilate
    return QC::Level::BadDoNotAssimilate;
  } else if (nDoNotAssimilate == length) {
    // No good observations, some mix of good and bad do-not-assimilate types
    //     -> GoodDoNotAssimilate
    return QC::Level::GoodDoNotAssimilate;
  } else {
    // No good observations, some mix of bad and do-not-assimilate -> Bad
    return QC::Level::Bad;
  }
}

void wholeReportFromPerProfile(const Data<feedback_io::QC::Level>& QCData,
    Data<feedback_io::QC::Level>& wholeReportQCData) {
  for (size_t iProfile = 0;
      iProfile < QCData.n_obs(); ++iProfile) {
    size_t nGoodObs = 0;
    size_t nBadObs = 0;
    size_t nBadDoNotAssimilate = 0;
    size_t nDoNotAssimilate = 0;
    for (size_t iLevel = 0;
        iLevel < QCData.length(iProfile);
        ++iLevel) {
      if (QCData(iProfile, iLevel) == QC::Level::Bad) {
        ++nBadObs;
      } else if (QCData(iProfile, iLevel) ==
                 QC::Level::BadDoNotAssimilate) {
        ++nBadDoNotAssimilate;
      }
      if (QCData(iProfile, iLevel) == QC::Level::BadDoNotAssimilate ||
          QCData(iProfile, iLevel) == QC::Level::GoodDoNotAssimilate ||
          QCData(iProfile, iLevel) == QC::Level::DoNotAssimilate) {
        ++nDoNotAssimilate;
      }
      if (QCData(iProfile, iLevel) == QC::Level::Good) {
        ++nGoodObs;
      }
    }
    wholeReportQCData[iProfile] = wholeReportUpdate(
        wholeReportQCData[iProfile], QCData.length(iProfile), nGoodObs,
        nBadObs, nBadDoNotAssimilate, nDoNotAssimilate);
  }
}
}  // namespace feedback_io
}  // namespace nemo_feedback
