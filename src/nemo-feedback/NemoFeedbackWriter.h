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
    size_t n_levels;
  };

  struct NameData {
    std::vector<std::string> variable_names;
    std::vector<std::string> additional_names;
    std::vector<std::string> long_names;
    std::vector<std::string> unit_names;
  };

  NemoFeedbackWriter(
      eckit::PathName& filename,
      const size_t & n_obs_to_write,
      const std::vector<bool> & to_write,
      const CoordData & coords,
      const NameData & name_data,
      const std::vector<bool> & extra_vars,
      const std::vector<std::string>& station_types,
      const std::vector<std::string>& station_ids);

  void write_variable_surf(
      const std::string & variable_name,
      const std::vector<double>& data);

  void write_variable_profile(
      const std::string & variable_name,
      const std::vector<double>& data,
      const std::vector<size_t>& record_starts,
      const std::vector<size_t>& record_counts);

  void write_variable_surf_qc(
      const std::string & variable_name,
      const std::vector<int32_t>& data);

  void write_variable_surf_qc(
      const std::string & variable_name,
      const std::vector<int32_t>& data,
      const size_t flag_index);

  void write_variable(
      const std::string & variable_name,
      const std::vector<double>& data) {
    write_variable_surf(variable_name, data);
  }

  void write_variable_level_qc(
      const std::string & variable_name,
      const std::vector<int32_t>& data,
      const std::vector<size_t>& record_starts,
      const std::vector<size_t>& record_counts);

  void write_variable_level_qc(
      const std::string & variable_name,
      const std::vector<int32_t>& data,
      const size_t flag_index,
      const std::vector<size_t>& record_starts,
      const std::vector<size_t>& record_counts);

  static constexpr double double_fillvalue = 99999.0;

 private:
  NemoFeedbackWriter() : ncFile(), nobs_dim(), nlevels_dim(), coords_(),
                         name_data_(), n_obs_(), n_obs_to_write_(), to_write_()
                       {}

  template <typename T>
  std::vector<T> reduce_data(
      const std::vector<T> & data_in);

  void define_coord_variables(
      const size_t n_obs_vars,
      const size_t n_add_entries,
      const size_t n_extra);

  void write_coord_variables();

  void write_metadata_variables(
      const std::vector<bool>& extra_vars);

  void define_whole_report_variables();

  void write_whole_report_variables(
      const std::vector<std::string> & station_types,
      const std::vector<std::string> & station_ids);

  void define_variable(
      const std::string & variable_name,
      const std::string & long_name,
      const std::string & unit_name);

  void define_extra_variable(
      const std::string & variable_name,
      const std::string & long_name,
      const std::string & unit_name);

  std::unique_ptr<netCDF::NcFile> ncFile;
  std::unique_ptr<netCDF::NcDim> nobs_dim;
  std::unique_ptr<netCDF::NcDim> nlevels_dim;
  std::unique_ptr<netCDF::NcDim> nqcf_dim;
  const CoordData coords_;
  const NameData name_data_;
  const size_t n_obs_;
  const size_t n_obs_to_write_;
  const std::vector<bool> to_write_;
  static const std::map<std::string, size_t> coord_sizes;
};
}  // namespace nemo_feedback

//------------------------------------------------------------------------------------------------------
