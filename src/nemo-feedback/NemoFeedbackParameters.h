/*
 * (C) Crown Copyright 2020, the Met Office. All rights reserved.
 *
 * Refer to COPYRIGHT.txt of this distribution for details.
 */

#ifndef OPSINPUTS_VAROBSWRITERPARAMETERS_H_
#define OPSINPUTS_VAROBSWRITERPARAMETERS_H_

#include <string>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "eckit/config/Configuration.h"
#include "oops/generic/ObsFilterParametersBase.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace nemo_feedback {

/// \brief NemoFeedback options.
class NemoFeedbackParameters : public oops::ObsFilterParametersBase {
  OOPS_CONCRETE_PARAMETERS(NemoFeedbackParameters, oops::ObsFilterParametersBase)

 public:
  oops::Parameter<std::string> Filename{"filename", "nemo_fdbk_out.nc", this};
  oops::Parameter<std::vector<eckit::LocalConfiguration>> Variables{"variables", {}, this};

};

}  // namespace nemo_feedback

#endif  // OPSINPUTS_VAROBSWRITERPARAMETERS_H_
