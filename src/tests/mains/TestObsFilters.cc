/*
 * (C) Crown Copyright 2023 Met Office
 */

#include "../test/ufo/ObsFilters.h"
#include "oops/runs/Run.h"
#include "nemo-feedback/instantiateObsFilterFactory.h"
#include "ufo/instantiateObsFilterFactory.h"
#include "ufo/ObsTraits.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  ufo::instantiateObsFilterFactory();
  nemo_feedback::instantiateObsFilterFactory<ufo::ObsTraits>();
  ufo::test::ObsFilters tests;
  return run.execute(tests);
}
