/*
 * (C) British Crown Copyright 2022 Met Office
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

#include "nemo-feedback/NemoFeedbackReduce.h"

namespace nemo_feedback {

/// \brief Interface to the NetCDF library to write feedback files
class NemoFeedbackWriter {
 public:
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
  };

  NemoFeedbackWriter(
      eckit::PathName& filename,
      const CoordData & coords,
      const NameData & name_data,
      const std::vector<bool> & extra_vars,
      const std::vector<std::string>& station_types,
      const std::vector<std::string>& station_ids);

  /// \brief Write surface variable data
  void write_variable_surf(
      const std::string & variable_name,
      const std::vector<double>& data);

  /// \brief Write profile variable data
  void write_variable_profile(
      const std::string & variable_name,
      const std::vector<double>& data);

  /// \brief Write surface QC data variable
  void write_variable_surf_qc(
      const std::string & variable_name,
      const std::vector<int32_t>& data);

  /// \brief Write surface QC data variable with specified flag
  void write_variable_surf_qc(
      const std::string & variable_name,
      const std::vector<int32_t>& data,
      const size_t flag_index);

  /// \brief Write level QC data variable
  void write_variable_level_qc(
      const std::string & variable_name,
      const std::vector<int32_t>& data);

  /// \brief Write level QC data variable with specified flag
  void write_variable_level_qc(
      const std::string & variable_name,
      const std::vector<int32_t>& data,
      const size_t flag_index);

  static constexpr double double_fillvalue = 99999.0;
  static constexpr int32_t int32_fillvalue = 0;

 private:
  NemoFeedbackWriter() : ncFile(), nobs_dim(), nlevels_dim(), coords_(),
                         name_data_() {}

  /// \brief Define the coordinate variables in the NetCDF file
  void define_coord_variables(
      const size_t n_obs_vars,
      const size_t n_add_entries,
      const size_t n_extra);

  /// \brief Write the coordinate variables to the NetCDF file
  void write_coord_variables();

  /// \brief Write the required feedback file metadata to the NetCDF file
  void write_metadata_variables(
      const std::vector<bool>& extra_vars);

  /// \brief Define the variables that impact entire observations in the
  ///        NetCDF file
  void define_whole_report_variables();

  /// \brief Write the variables that impact entire observations in the
  ///        NetCDF file
  void write_whole_report_variables(
      const std::vector<std::string> & station_types,
      const std::vector<std::string> & station_ids);

  /// \brief Define a data variable in the NetCDF file
  void define_variable(
      const std::string & variable_name,
      const std::string & long_name,
      const std::string & unit_name);

  /// \brief Define an 'extra' data variable in the NetCDF file
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
  static const std::map<std::string, size_t> coord_sizes;
};
}  // namespace nemo_feedback

//------------------------------------------------------------------------------------------------------
