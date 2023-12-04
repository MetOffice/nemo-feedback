/*
 * (C) British Crown Copyright 2023 Met Office
 */

#include "nemo-feedback/feedback_io/Utils.h"

#include <netcdf>
// https://github.com/Unidata/netcdf-cxx4

#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cstring>

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"
#include "oops/util/Duration.h"

#define DOUBLE_FILLVALUE 99999.0
#define FLOAT_FILLVALUE 99999.0
#define INT32_FILLVALUE 0
#define STRING_FILLVALUE "MISSING "


namespace nemo_feedback {
namespace feedback_io {

// default
template <typename T> const netCDF::NcType NetCDFTypeMap<T>::ncType =
  netCDF::ncDouble;

template <> const netCDF::NcType NetCDFTypeMap<float>::ncType = netCDF::ncFloat;

template const netCDF::NcType NetCDFTypeMap<float>::ncType;
template const netCDF::NcType NetCDFTypeMap<double>::ncType;

template <> const double typeToFill::value<double>() {
  return DOUBLE_FILLVALUE;
}
template <> const float  typeToFill::value<float>() {
  return FLOAT_FILLVALUE;
}
template <> const int32_t  typeToFill::value<int32_t>() {
  return INT32_FILLVALUE;
}
template <> const QC::Level  typeToFill::value<QC::Level>() {
  return QC::Level::None;
}
template <> const std::string  typeToFill::value<std::string>() {
  return STRING_FILLVALUE;
}

void NameData::validate() const {
  const size_t n_vars = variable_names.size();
  std::string className = "nemo_feedback::feedback_io::NameData";
  ASSERT_MSG(legacy_ops_qc_conventions.size() == n_vars,
      className + " legacy_ops_qc_conventions.size() mismatch "
      + "with number of variables ");
  ASSERT_MSG(long_names.size() == n_vars,
      className + "long_names.size() size mismatch with number of variables ");
  ASSERT_MSG(unit_names.size() == n_vars,
      className + "unit_names.size() mismatch with number of variables ");
}
}  // namespace feedback_io
}  // namespace nemo_feedback
