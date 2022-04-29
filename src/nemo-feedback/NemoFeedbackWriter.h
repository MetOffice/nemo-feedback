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

  struct CoordData {
    std::vector<double> lats;
    std::vector<double> lons;
    std::vector<double> depths;
    std::vector<double> julian_days;
    util::DateTime juld_reference;
  };

  NemoFeedbackWriter(
      eckit::PathName& filename,
      const size_t & n_obs_to_write,
      const std::vector<bool> & to_write,
      const std::vector<double> & lons,
      const std::vector<double> & lats,
      const std::vector<double> & depths,
      const std::vector<double> & times,
      const std::vector<std::string> & variable_names,
      const std::vector<std::string> & long_names,
      const std::vector<std::string> & unit_names,
      const std::vector<std::string> & additional_variables,
      const std::vector<bool> & extra_vars,
      const size_t n_levels,
      const util::DateTime & juld_reference,
      const std::vector<std::string>& station_types,
      const std::vector<std::string>& station_ids);

  void write_variable_surf(
      const size_t & n_obs,
      const size_t & n_obs_to_write,
      const std::vector<bool> & to_write,
      const std::string & variable_name,
      const std::vector<double>& data);

  void write_variable_profile(
      const size_t & n_obs,
      const std::string & variable_name,
      const std::vector<double>& data,
      const std::vector<size_t>& record_starts,
      const std::vector<size_t>& record_counts);

  void write_variable_surf_qc(
      const size_t & n_obs,
      const size_t & n_obs_to_write,
      const std::vector<bool> & to_write,
      const std::string & variable_name,
      const std::vector<int32_t>& data);

  void write_variable_surf_qc(
      const size_t & n_obs,
      const size_t & n_obs_to_write,
      const std::vector<bool> & to_write,
      const std::string & variable_name,
      const std::vector<int32_t>& data,
      const size_t flag_index);

  void write_variable(
      const size_t & n_obs,
      const size_t & n_obs_to_write,
      const std::vector<bool> & to_write,
      const std::string & variable_name,
      const std::array<double, 2>& data) {}

  void write_variable_level_qc(
      const size_t & n_obs,
      const size_t & n_obs_to_write,
      const std::vector<bool> & to_write,
      const std::string & variable_name,
      const std::array<int32_t, 2>& data) {}

  static constexpr double double_fillvalue = 99999.0;

 private:
  NemoFeedbackWriter() : ncFile(), nobs_dim(), nlevels_dim(), n_levels_() {}

  template <typename T>
  std::vector<T> reduce_data(
      const size_t & n_obs,
      const size_t & n_obs_to_write,
      const std::vector<bool> & to_write,
      const std::vector<T> & data_in);

  void define_coord_variables(
      const size_t n_obs_to_write,
      const size_t n_levels,
      const size_t n_obs_vars,
      const size_t n_add_entries,
      const size_t n_extra);

  void write_coord_variables(
      const size_t & n_obs,
      const size_t & n_obs_to_write,
      const std::vector<bool> & to_write,
      const std::vector<double>& lons,
      const std::vector<double>& lats,
      const std::vector<double>& levels,
      const std::vector<double>& times);

  void write_metadata_variables(
      const std::vector<std::string>& variable_names,
      const std::vector<std::string> & additional_variables,
      const std::vector<bool>& extra_vars,
      const util::DateTime& juld_reference);

  void define_whole_report_variables();

  void write_whole_report_variables(
      const size_t & n_obs,
      const size_t & n_obs_to_write,
      const std::vector<bool> & to_write,
      const std::vector<std::string> & station_types,
      const std::vector<std::string> & station_ids);

  void define_variable(
      const std::string & variable_name,
      const std::string & long_names,
      const std::string & unit_names,
      const std::vector<std::string> & additional_names =
        std::vector<std::string>());

  void define_extra_variable(
      const std::string & variable_name,
      const std::string & long_name,
      const std::string & unit_name);

  std::unique_ptr<netCDF::NcFile> ncFile;
  std::unique_ptr<netCDF::NcDim> nobs_dim;
  std::unique_ptr<netCDF::NcDim> nlevels_dim;
  std::unique_ptr<netCDF::NcDim> nqcf_dim;
  size_t n_levels_;
  static const std::map<std::string, size_t> coord_sizes;
};
}  // namespace nemo_feedback

//------------------------------------------------------------------------------------------------------
