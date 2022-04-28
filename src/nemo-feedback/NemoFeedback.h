/*
 * (C) Crown Copyright 2021, the Met Office. All rights reserved.
 *
 * Refer to COPYRIGHT.txt of this distribution for details.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "ioda/ObsDataVector.h"
#include "oops/base/Variables.h"
#include "oops/interface/ObsFilterBase.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "nemo-feedback/NemoFeedbackParameters.h"
#include "ufo/filters/ObsFilterData.h"
#include "ufo/ObsTraits.h"

namespace nemo_feedback {

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
 void setupJulianDay(const size_t n_obs,
                     util::DateTime& juld_reference,
                     std::vector<double>& julian_days) const;
  void setupAltimeterIds(const size_t n_obs,
                        std::vector<std::string>& station_ids,
                        std::vector<std::string>& station_types,
                        std::vector<bool>& to_write) const;
  void setupSurfaceIds(const size_t n_obs,
                       std::vector<std::string>& station_ids,
                       std::vector<std::string>& station_types) const;
  void print(std::ostream &) const override;

  ioda::ObsSpace & obsdb_;
  ufo::ObsFilterData data_;
  oops::Variables geovars_;
  oops::Variables extradiagvars_;
  std::shared_ptr<ioda::ObsDataVector<int>> flags_;
  std::shared_ptr<ioda::ObsDataVector<float>> obsErrors_;
  const util::DateTime validityTime_;

  NemoFeedbackParameters parameters_;
};

}  // namespace nemo_feedback

