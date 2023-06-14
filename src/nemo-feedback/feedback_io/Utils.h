/*
 * (C) British Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#pragma once

#include <netcdf>

#include <string>
#include <memory>
#include <vector>

#include "oops/util/DateTime.h"

#define DOUBLE_FILLVALUE 99999.0
#define FLOAT_FILLVALUE 99999.0
#define INT32_FILLVALUE 0
#define STRING_FILLVALUE "MISSING "

namespace nemo_feedback {
namespace feedback_io {

// create a mapping between C++ types and NetCDF type objects
template <typename T>
struct NetCDFTypeMap {
  static const netCDF::NcType ncType;
};

// create a mapping between C++ types and Nemo fillvalues
class typeToFill {
 public:
  template <typename T> static const T value();
};

/// \brief Coordinate information for the feedback file, vectors with length
///        N_OBS, apart from depths which must have dimension N_OBS*N_LEVELS
struct CoordData {
  std::vector<double> lats;
  std::vector<double> lons;
  std::vector<double> julian_days;
  std::vector<double> depths;
  std::vector<size_t> record_starts;
  std::vector<size_t> record_counts;
  util::DateTime juld_reference;
  size_t n_levels;
  size_t n_obs;
  size_t n_locs;
};

/// \brief Naming information for each variable in the feedback file
struct NameData {
  std::vector<std::string> variable_names;
  std::vector<std::string> additional_names;
  std::vector<std::string> long_names;
  std::vector<std::string> unit_names;
  std::vector<bool> legacy_ops_qc_conventions;

  void validate() const;
};
}  // namespace feedback_io
}  // namespace nemo_feedback
