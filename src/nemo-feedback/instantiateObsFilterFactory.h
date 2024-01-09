/*
 * (C) British Crown Copyright 2024 Met Office
 */

#pragma once

#include "oops/interface/ObsFilterBase.h"
#include "nemo-feedback/NemoFeedback.h"

namespace nemo_feedback {

template<typename OBS>
void instantiateObsFilterFactory() {
  static oops::interface::FilterMaker<OBS, NemoFeedback>
    makerNemoFeedback_("NEMO Feedback Writer");
}

}  // namespace nemo_feedback

