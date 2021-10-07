/*
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#pragma once

#include <netcdf>

#include <string>
#include <memory>
#include <map>
#include <vector>

#include "eckit/filesystem/PathName.h"

#include "oops/util/DateTime.h"

#include "atlas/runtime/Exception.h"

namespace nemo_feedback {

class NemoFeedbackWriter {
 public:
  NemoFeedbackWriter(eckit::PathName& filename,
      const std::vector<double>& lons, const std::vector<double>& lats,
      const std::vector<double>& depths, const std::vector<double>& times,
      const std::vector<std::string>& variable_names,
      const std::vector<std::string>& additional_variables,
      const size_t n_levels, const util::DateTime & juld_reference);

  void write_variable_surf(const std::string & variable_name,
      const std::vector<double>& data);

  void write_variable_surf_qc(const std::string & variable_name,
      const std::vector<int32_t>& data);
  void write_variable_surf_qc(const std::string & variable_name,
      const std::vector<int32_t>& data, const int32_t flag_index);

  void write_variable(const std::string & variable_name,
      const std::array<double, 2>& data) {}

  void write_variable_level_qc(const std::string & variable_name,
      const std::array<int32_t, 2>& data) {}

  static constexpr double double_fillvalue = 99999.0;

 private:
  NemoFeedbackWriter() : ncFile(), nobs_dim(), nlevels_dim(), n_levels_() {}

  void define_coord_variables(const size_t n_obs, const size_t n_levels,
      const size_t n_obs_vars, const size_t n_add_entries);

  void write_coord_variables(const std::vector<double>& lons,
      const std::vector<double>& lats, const std::vector<double>& levels,
      const std::vector<double>& times);

  void write_metadata_variables(
      const std::vector<std::string>& variable_names,
      const std::vector<std::string> & additional_variables,
      const util::DateTime& juld_reference);

  void define_whole_report_variables(const size_t n_obs,
      const size_t n_levels);

  void define_variable(const std::string & variable_name,
      const std::vector<std::string> & additional_names =
        std::vector<std::string>());

  std::unique_ptr<netCDF::NcFile> ncFile;
  std::unique_ptr<netCDF::NcDim> nobs_dim;
  std::unique_ptr<netCDF::NcDim> nlevels_dim;
  std::unique_ptr<netCDF::NcDim> nqcf_dim;
  size_t n_levels_;
  static const std::map<std::string, size_t> coord_sizes;
};
}  // namespace nemo_feedback

//------------------------------------------------------------------------------------------------------
