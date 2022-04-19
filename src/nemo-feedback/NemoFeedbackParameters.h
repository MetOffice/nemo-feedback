/*
 * (C) Crown Copyright 2020, the Met Office. All rights reserved.
 *
 * Refer to COPYRIGHT.txt of this distribution for details.
 */

#pragma once

#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "oops/generic/ObsFilterParametersBase.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/RequiredParameter.h"

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
  /// Logic used to select locations to be written to file.
  /// If not specified, all locations will be written.
  oops::Parameter<std::vector<ufo::WhereParameters>> where{"where", {}, this};
};

}  // namespace nemo_feedback

