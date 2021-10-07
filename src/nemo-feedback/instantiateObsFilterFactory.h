/*
 * (C) Crown Copyright 2020, the Met Office. All rights reserved.
 *
 * Refer to COPYRIGHT.txt of this distribution for details.
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

