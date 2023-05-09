/*
 * (C) Crown Copyright 2021, the Met Office. All rights reserved.
 *
 * Refer to COPYRIGHT.txt of this distribution for details.
 */

#pragma once
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/mpi/Comm.h"
#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/interface/ObsFilterBase.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "nemo-feedback/NemoFeedbackParameters.h"
#include "nemo-feedback/NemoFeedbackWriter.h"
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
  void write_all_data(NemoFeedbackWriter<T>& fdbk_writer,
                      NemoFeedbackReduce& reducer,
                      const CoordData& coords,
                      const ioda::ObsVector &ov,
                      const bool is_profile) const;
  /// \brief group coordinate data by the ioda record whilst filtering
  ///        observations to write to the file
  void groupCoordsByRecord(const std::vector<bool>& to_write,
                          CoordData& coords,
                          bool is_profile) const;
  /// \brief Setup the NEMO STATION_TYPES and STATION_IDS netCDF variables for
  ///        altimeter observations.  Filter observations based on the latest
  ///        version. TODO: move this into another UFO filter.
  void setupAltimeterIds(const size_t n_obs,
                        std::vector<std::string>& station_ids,
                        std::vector<std::string>& station_types,
                        std::vector<bool>& to_write) const;
  /// \brief Setup the NEMO STATION_TYPES and STATION_IDS netCDF variables
  void setupIds(const size_t n_obs,
                const std::vector<size_t>& record_starts,
                const std::vector<size_t>& record_counts,
                std::vector<std::string>& station_ids,
                std::vector<std::string>& station_types) const;
  /// \brief Sync relevent coordinate data across MPI processes
  void mpi_sync_coordinates(CoordData& coords,
                            const eckit::mpi::Comm& comm);
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
};

}  // namespace nemo_feedback

