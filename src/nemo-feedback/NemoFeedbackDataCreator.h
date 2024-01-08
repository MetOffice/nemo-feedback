/*
 * (C) British Crown Copyright 2024 Met Office
 */

#pragma once
#include <memory>
#include <string>
#include <vector>
#include <tuple>

#include "eckit/mpi/Comm.h"
#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "ufo/filters/DiagnosticFlag.h"
#include "nemo-feedback/feedback_io/Data.h"
#include "nemo-feedback/feedback_io/DataIndexer.h"

namespace nemo_feedback {

/// \brief Convert ioda/ufo data to FeedbackData objects
class NemoFeedbackDataCreator {
 public:
  static const std::string className()
    {return "nemo_feedback::NemoFeedbackDataCreator";}

  /// \brief setup the indexer on construction using ioda profile information
  NemoFeedbackDataCreator(const ioda::ObsSpace& obsdb,
      const ioda::ObsVector& hofxObsVector);

  /// \brief setup the indexer on construction using ioda profile information
  ///        and sub-select observations
  NemoFeedbackDataCreator(const ioda::ObsSpace& obsdb,
      const ioda::ObsVector& hofxObsVector, std::vector<bool>& select) :
    NemoFeedbackDataCreator(obsdb, hofxObsVector) {
    update_indexer(select);
  }

  /// \brief update indexer with new, reduced size indexer of a subset of the
  ///        data. Any objects created after this update will have the new
  ///        indexer.
  void update_indexer(std::vector<bool>& select) {
    indexer_ = std::make_shared<feedback_io::DataIndexer>((*indexer_),  select);
  }

  /// \brief create time feedbackData
  std::tuple<std::string, feedback_io::Data<double>> create_datetimes(
     const std::string& obsGroup, const std::string& ufoName,
     const util::DateTime juldReferenceDT) const;

  /// \brief create string feedbackData (station IDs and types)
  template<typename T>
  feedback_io::Data<std::string> create(const std::string& obsGroup,
      const std::string& ufoName,
      const T typeInstance,
      size_t width, bool leftJustify = false) const;

  /// \brief create altimiter id feedback data
  std::tuple<feedback_io::Data<std::string>, feedback_io::Data<std::string>>
    create_altimeter_IDs() const;

  /// \brief create DiagnosticFlags feedbackData
  feedback_io::Data<feedback_io::QC::Level> create(const std::string& obsGroup,
      const std::string& ufoName,
      const ufo::DiagnosticFlag typeInstance,
      const feedback_io::QC::Level whenTrue,
      const feedback_io::QC::Level whenFalse) const;

  /// \brief create feedbackData
  template<typename T>
  feedback_io::Data<T> create(const std::string& obsGroup,
      const std::string& ufoName,
      const T typeInstance) const;

std::shared_ptr<feedback_io::DataIndexer> indexer() const {return indexer_;}

 private:
  /// \brief create station type feedbackData from obsdb
  feedback_io::Data<std::string> create_from_obsdb(const std::string& obsGroup,
      const std::string& ufoName,
      const int32_t typeInstance,
      size_t width, bool leftJustify = false) const;

  /// \brief create station ID feedbackData from obsdb
  feedback_io::Data<std::string> create_from_obsdb(const std::string& obsGroup,
      const std::string& ufoName,
      const std::string typeInstance,
      size_t width, bool leftJustify = false) const;

  /// \brief create DiagnosticFlags feedbackData from obsdb
  feedback_io::Data<feedback_io::QC::Level> create_from_obsdb(
      const std::string& obsGroup,
      const std::string& ufoName,
      const ufo::DiagnosticFlag TypeInstance,
      feedback_io::QC::Level whenTrue,
      feedback_io::QC::Level whenFalse) const;

  /// \brief create feedbackData from obsdb
  template<typename T>
  feedback_io::Data<T> create_from_obsdb(const std::string& obsGroup,
      const std::string& ufoName,
      const T typeInstance) const;

  /// \brief create feedbackData from hofx
  template<typename T>
  feedback_io::Data<T> create_from_hofx(const std::string& ufoName,
      const T typeInstance) const;

  const ioda::ObsSpace& obsdb_;
  const ioda::ObsVector& hofxObsVector_;
  std::shared_ptr<feedback_io::DataIndexer> indexer_;
};

}  // namespace nemo_feedback

