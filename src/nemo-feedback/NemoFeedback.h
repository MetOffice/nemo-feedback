/*
 * (C) British Crown Copyright 2024 Met Office
 */

#pragma once
#include <memory>
#include <ostream>
#include <string>
#include <vector>
#include <tuple>

#include "eckit/mpi/Comm.h"
#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/interface/ObsFilterBase.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "nemo-feedback/NemoFeedbackParameters.h"
#include "nemo-feedback/feedback_io/Writer.h"
#include "nemo-feedback/NemoFeedbackDataCreator.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/ObsTraits.h"
#include "ufo/utils/VariableNameMap.h"

namespace nemo_feedback {

/// \brief UFO filter for outputting data to NEMO feedback file
class NemoFeedback : public oops::interface::ObsFilterBase<ufo::ObsTraits>,
                     private util::ObjectCounter<NemoFeedback> {
 public:
  static const std::string classname() {return "nemo_feedback::NemoFeedback";}

  /// The type of parameters accepted by the constructor of this filter.
  /// This typedef is used by the FilterFactory.
  typedef NemoFeedbackParameters Parameters_;

  NemoFeedback(ioda::ObsSpace &, const Parameters_ &,
               std::shared_ptr<ioda::ObsDataVector<int> > flags,
               std::shared_ptr<ioda::ObsDataVector<float> > obsErrors);
  ~NemoFeedback();

  void preProcess() override {}
  void priorFilter(const ufo::GeoVaLs &) override;
  void postFilter(const ufo::GeoVaLs & gv,
                  const ioda::ObsVector &ov,
                  const ioda::ObsVector &bv,
                  const ufo::ObsDiagnostics &dv) override;
  void checkFilterData(const oops::FilterStage filterStage) override {}

  oops::Variables requiredVars() const override {return geovars_;}
  oops::Variables requiredHdiagnostics() const override {return extradiagvars_;}

 private:
  /// \brief write the data to the feedback file depending on chosen type
  template <typename T>
  void write_all_data(feedback_io::Writer<T>& fdbk_writer,
                      const NemoFeedbackDataCreator& creator) const;
  /// \brief Update the observations to write based upon the latest available
  ///        version of the data. This step ought to be moved to a separate ufo
  ///        filter when possible.
  void updateAltimeterSelection(std::vector<bool>& to_write) const;
  /// \brief setup the feedback file metadata - this is data that is written to
  ///        the file during construction of the feedback file. This step should
  ///        ensure that the number of levels and the reference Julian day are
  ///        shared across MPI processes.
  feedback_io::MetaData setupMetaData(const NemoFeedbackDataCreator& creator)
    const;
  /// \brief Setup the NEMO STATION_TYPES and STATION_IDS netCDF variables
  ///        In the case of altimetry data, this requires additional processing
  ///        of the satelliteIdentifier to extract these variables.
  std::tuple<feedback_io::Data<std::string>, feedback_io::Data<std::string>>
    setupIDs(const NemoFeedbackDataCreator& creator) const;
  /// \brief Sync number of levels and reference julian day data across MPI
  /// processes
  std::tuple<util::DateTime, size_t> mpiSync(size_t nLevelsLocal) const;
  void print(std::ostream &) const override;

  ioda::ObsSpace & obsdb_;
  ufo::ObsFilterData data_;
  oops::Variables geovars_;
  oops::Variables extradiagvars_;
  std::shared_ptr<ioda::ObsDataVector<int>> flags_;
  std::shared_ptr<ioda::ObsDataVector<float>> obsErrors_;
  NemoFeedbackParameters parameters_;
  ufo::VariableNameMap nameMap_;
  const util::DateTime validityTime_;
  feedback_io::NameData nameData_;
  std::vector<bool> isExtraVariable_;
  bool isAltimeter_;
};

}  // namespace nemo_feedback

