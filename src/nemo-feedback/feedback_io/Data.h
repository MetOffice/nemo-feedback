/*
 * (C) British Crown Copyright 2023 Met Office
 */


#pragma once

#include <nemo-feedback/feedback_io/DataIndexer.h>

#include <string>
#include <vector>
#include <memory>

#include <nemo-feedback/feedback_io/Utils.h>

namespace nemo_feedback {
namespace feedback_io {

  /// \brief Package a 1D vector of data with its indexer to the sorted
  ///        (profile, level) or (location) for surface information
template<typename T>
class Data {
 public:
  static std::string className() { return "nemo_feedback::feedback_io::Data"; }

  /// \brief empty constructor
  Data() :
    indexer_(std::make_shared<DataIndexer>()),
    data_(std::make_shared<std::vector<T>>(std::vector<T>())) {}

  /// \brief basic new constructor with a shallow copy of the data
  ///        (only wraps pointers to supplied objects)
  Data(const std::shared_ptr<DataIndexer>& indexer,
                   const std::shared_ptr<std::vector<T>> data) :
    indexer_(indexer), data_(data) { validate(); }

  /// \brief basic new constructor with a deep copy of the data
  ///        (creates new pointer with deep copy of the data)
  Data(const std::shared_ptr<DataIndexer>& indexer,
                   const std::vector<T>& data) :
    indexer_(indexer), data_(std::make_shared<std::vector<T>>(data))
    { validate(); }

  /// \brief basic new constructor with a deep copy of the data and indexer
  ///        (creates new pointer with deep copy of the data and indexer)
  Data(const DataIndexer& indexer,
                   const std::vector<T>& data) :
    indexer_(std::make_shared<DataIndexer>(indexer)),
    data_(std::make_shared<std::vector<T>>(data)) { validate(); }

  /// \brief default copy constructor
  Data(const Data&) = default;

  /// \brief reduced shallow copy
  ///   this might not be the most efficient as creates a new indexer for every
  ///   data object, prefer update indexer where possible
  template <typename B>
  Data(const Data<B>& other, const std::vector<bool>& select) :
    indexer_(std::make_shared<DataIndexer>((*other.indexer_), select)),
    data_(other.data_) { validate(); }

  /// \brief deep copy constructor
  ///   Compresses the source data whilst copying
  Data<T> deep_copy() const;

  /// \brief update indexer to point to indexer in other data object
  template <typename B>
  void update_indexer(const Data<B>& other) {
    indexer_ = other.indexer();
    validate();
  }

  /// \brief update indexer to point to other indexer pointer
  void update_indexer(const std::shared_ptr<DataIndexer> indexer) {
    indexer_ = indexer;
    validate();
  }

  /// \brief update indexer to new reduced indexer
  void update_indexer(const std::vector<bool>& select) {
    indexer_ = std::make_shared<DataIndexer>((*indexer_), select);
    validate();
  }

  /// \brief index into surface data indexed via nObs dimension
  T& operator[](size_t iLocation) { return (*data_)[(*indexer_)(iLocation)]; }
  T operator[](size_t iLocation) const
    { return (*data_)[(*indexer_)(iLocation)]; }
  /// \brief index into profile data indexed via nObs, nLevels dimensions
  T& operator()(size_t iProfile, size_t iLevel)
    { return (*data_)[(*indexer_)(iProfile, iLevel)]; }
  T operator()(size_t iProfile, size_t iLevel) const
    { return (*data_)[(*indexer_)(iProfile, iLevel)]; }
  /// \brief index into surface data with bounds-checking indexed via nObs
  ///        dimension
  T at(size_t iLocation) const { return data_->at(indexer_->at(iLocation)); }
  /// \brief index into profile data with bounds-checking indexed via nObs,
  ///        nLevels dimensions
  T at(size_t iProfile, size_t iLevel) const
    { return data_->at(indexer_->at(iProfile, iLevel)); }
  /// \brief retrieve all valid data
  std::vector<T> raw() const;
  /// \brief retrieve all valid data at surface
  std::vector<T> raw_surface() const;
  /// \brief retrieve all valid data for a profile
  std::vector<T> raw_profile(size_t iProfile) const;
  /// \brief number of elements in a profile
  size_t length(size_t iProfile) const { return indexer_->length(iProfile); }
  /// \brief number of observations in data
  size_t n_obs() const { return indexer_->n_obs(); }
  /// \brief maximum number of levels across data
  size_t n_levels() const { return indexer_->n_levels(); }
  size_t n_locations() const { return indexer_->n_locations(); }
  const std::shared_ptr<DataIndexer> indexer() const { return indexer_; }

 private:
  void validate() const;
  std::shared_ptr<DataIndexer> indexer_;
  std::shared_ptr<std::vector<T>> data_;
};

/// \brief create whole report data from per-profile QC data
void wholeReportFromPerProfile(const Data<feedback_io::QC::Level>& QCData,
    Data<feedback_io::QC::Level>& wholeReportQCData);
}  // namespace feedback_io
}  // namespace nemo_feedback
