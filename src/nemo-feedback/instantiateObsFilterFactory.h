/*
 * (C) British Crown Copyright 2024 Met Office
 */

#pragma once

#include "nemo-feedback/NemoFeedback.h"
#include "ufo/ObsFilterBase.h"

namespace nemo_feedback {

void instantiateObsFilterFactory() {
  static ufo::FilterMaker<NemoFeedback>
    makerNemoFeedback_("NEMO Feedback Writer");
}

}  // namespace nemo_feedback

