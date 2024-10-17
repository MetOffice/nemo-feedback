/*
 * (C) British Crown Copyright 2024 Met Office
 */

#include "../test/ufo/ObsFilters.h"
#include "oops/runs/Run.h"
#include "nemo-feedback/instantiateObsFilterFactory.h"
#include "ufo/instantiateObsFilterFactory.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  ufo::instantiateObsFilterFactory();
  nemo_feedback::instantiateObsFilterFactory();
  ufo::test::ObsFilters tests;
  return run.execute(tests);
}
