![CI](https://github.com/MetOffice/nemo-feedback/actions/workflows/ci.yml/badge.svg)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
&copy; British Crown Copyright 2023, the Met Office. All rights reserved.

# nemo-feedback

JEDI UFO filter extension for writing NEMO ocean model NetCDF feedback file format data from the UK Met Office.

## Description

_nemo-feedback_ is a small code snippet that is compiled with UFO to add an additional filter for direct output in NEMO feedback file format. It was developed at the Met Office to interface with the JEDI framework and the [_orca-jedi_](https://github.com/MetOffice/nemo-feedback) JEDI "psuedomodel". These tools are built as a bridge between NEMO/NEMOVAR and the JEDI framework - rather than interfacing directly with NEMO or NEMOVAR, the model state is derived from input files, and the quality controlled observations are written to feedback file for ingestion by the NEMOVAR executable. It is anticipated that _nemo-feedback_ will become obsolete, either through adoption of a standard observation file format across Met Office systems (ODB), or through further integration of Met Office Ocean workflows into the JEDI framework.

## Getting Started

### Dependencies

  * [cmake](https://cmake.org/)
  * [Unidata/netcdf-cxx4](https://github.com/Unidata/netcdf-cxx4)
  * [ecmwf/ecbuild](https://github.com/ecmwf/ecbuild)
  * [ecmwf/eckit](https://github.com/ecmwf/eckit)
  * [JCSDA/oops](https://github.com/JCSDA/oops)
  * [JCSDA/ufo](https://github.com/JCSDA/ufo)

### Installing

_nemo-feedback_ meant to be built with UFO as part of [mo-bundle](https://github.com/MetOffice/mo-bundle) at the Met Office - see the README in that project for details on how to build.

Otherwise, it is possible to build _nemo-feedback_ within a custom JEDI bundle. For details about JEDI, including installation instructions see the [jedi-docs](http://jedi-docs.jcsda.org/). A small example of this is part of the continuous integration in this repository.

These bundles are built, made and installed via cmake, and tested with ctest. All code should be documented at the source level for processing using doxygen.

### Using the interface

The jedi configuration is documented in the main jedi documentation. Settings for Met Office operational numerical weather prediction workflows are held internally by the Met Office, however there are some example configurations in the ``examples`` directory as well as the ctest tests inputs (``src/tests/testinputs``). The _nemo-feedback_ filter will work both serially or when using MPI.

The parameters for _nemo-feedback_ are documented in code in the parameters class. When the program is compiled, these are exported to a yaml schema json file, that can be used in conjunction with your editor to highlight any issues with your configuration. An example yaml configuration would look something like:

```yaml
...
observations:
  ...
  - obs space:
    ...
    obs filters:
    ...
    - filter: NEMO Feedback Writer
      reference date:  # eckit::DateTime reference date for the julian day observation time in the file
      filename: # path to output file
      observation alias file:  # mapping between observation and model names inside UFO
      variables:  # List of variables to write
      - name:  # UFO observation variable name
        nemo name:  # feedback file variable name
        long name:  # long name in the netcdf file
        units: # units attribute in the netcdf file
        extra variable:  # optional, if true don't make the normal set of sub-variables (e.g for mean-dynamic-topography)
        additional variables: # "additional" variables (these are extra sub-variables of a variable)
        - name: # UFO observation variable name
          feedback suffix: # suffix for this additional variable name (i.e <nemo name>_<feedback suffix>)
          ioda group: # IODA group to extract this variable from
...
```

## Help

See the JEDI documentation for help. Additional debugging/trace output is available when:
```
OOPS_DEBUG=true
OOPS_TRACE=true
```

## Authors

The current lead maintainer is [@twsearle](https://github.com/twsearle) along with a large amount of help from Met Office contributors (see the "Contributors" page on github).

## Working practices

Please see the [JEDI working principles](https://jointcenterforsatellitedataassimilation-jedi-docs.readthedocs-hosted.com/en/latest/working-practices/index.html) for current working practises. There are also templates for issues and PRs should you wish to contribute.

