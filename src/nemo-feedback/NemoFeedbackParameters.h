/*
 * (C) British Crown Copyright 2024 Met Office
 */

#pragma once

#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "oops/generic/ObsFilterParametersBase.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "ufo/filters/processWhere.h"
#include "nemo-feedback/NemoFeedbackParameterTraitsOutputDtype.h"

namespace nemo_feedback {

/// \brief NemoFeedback options for an "additional" variable.
class NemoFeedbackAddVariableParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(NemoFeedbackAddVariableParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::string> name{"name", this};
  oops::RequiredParameter<std::string> feedbackSuffix{"feedback suffix", this};
  oops::RequiredParameter<std::string> iodaGroup{"ioda group", this};
};

/// \brief NemoFeedback options for a "main" variable.
class NemoFeedbackVariableParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(NemoFeedbackVariableParameters, oops::Parameters)

 public:
  oops::RequiredParameter<std::string> name{"name", this};
  oops::RequiredParameter<std::string> nemoName{"nemo name", this};
  oops::OptionalParameter<std::string> iodaObsGroup{"ioda group", this};
  oops::RequiredParameter<std::string> units{"units", this};
  oops::RequiredParameter<std::string> longName{"long name", this};
  oops::OptionalParameter<bool> extravar{"extra variable", this};
  oops::Parameter<std::vector<NemoFeedbackAddVariableParameters>>
    variables{"additional variables", {}, this};
};

/// \brief NemoFeedback options.
class NemoFeedbackParameters : public oops::ObsFilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(NemoFeedbackParameters,
      oops::ObsFilterParametersBase)

 public:
  oops::Parameter<std::string> Filename{"filename", "nemo_fdbk_out.nc", this};
  oops::RequiredParameter<std::vector<NemoFeedbackVariableParameters>>
    variables{"variables", this};
  oops::OptionalParameter<util::DateTime> refDate{"reference date", this};
  oops::OptionalParameter<std::string> depthGroup{"depth group", this};
  oops::OptionalParameter<std::string> depthVariable{"depth variable", this};
  /// Data Type (float or default to double) of the netCDF output data
  oops::OptionalParameter<OutputDtype> type{"type", this};
  /// \brief boolean mask to select locations to be written to file.
  ///        If not specified, all locations will be written.
  oops::Parameter<std::vector<ufo::WhereParameters>> where{"where", {}, this};
  /// \brief Parameter specifying path to yaml file containing Observation to
  ///        GeoVaL name mapping
  oops::OptionalParameter<std::string> geoVaLsAliasFile{
    "observation alias file", this};
};

}  // namespace nemo_feedback
