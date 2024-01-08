/*
 * (C) British Crown Copyright 2024 Met Office
 */


#pragma once

#include <netcdf>

#include <string>
#include <memory>
#include <map>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "eckit/filesystem/PathName.h"

#include "nemo-feedback/feedback_io/Utils.h"
#include "nemo-feedback/feedback_io/Data.h"

namespace nemo_feedback {
namespace feedback_io {

/// \brief MetaData for the feedback file containing the
/// coordinates and station information for each observation
struct MetaData {
 public:
  static std::string className()
  { return "nemo_feedback::feedback_io::MetaData"; }
  /// \note this constructor should never be used - it is here so that the
  /// explicit Writer() constructor can be made private and thus force
  /// initialisation at allocation. TODO: do the same for this constructor.
  MetaData() : nObs(0), nLevels(0), nLocations(0) {}
  MetaData(const Data<double>& latsIn,
           const Data<double>& lonsIn,
           const Data<double>& julianDaysIn,
           const Data<double>& depthsIn,
           const Data<std::string>& stationTypesIn,
           const Data<std::string>& stationIDsIn,
           size_t nLevelsGlobal,
           const std::string& juldReferenceIn) : lats(latsIn), lons(lonsIn),
    julianDays(julianDaysIn), depths(depthsIn), stationTypes(stationTypesIn),
    stationIDs(stationIDsIn), juldReference(juldReferenceIn),
    nObs(latsIn.n_obs()), nLevels(nLevelsGlobal),
    nLocations(latsIn.n_locations()) {
    this->validate();
  }
  void validate() {
    ASSERT_MSG(lons.n_locations() == nLocations,
        MetaData::className() + ": number of longitudes does not match");
    ASSERT_MSG(julianDays.n_locations() == nLocations,
        MetaData::className() + ": number of julian days does not match");
    ASSERT_MSG(depths.n_locations() == nLocations,
        MetaData::className() + ": number of depths does not match");
    ASSERT_MSG(stationTypes.n_locations() == nLocations,
        MetaData::className() + ": number of station types does not match");
    ASSERT_MSG(stationIDs.n_locations() == nLocations,
        MetaData::className() + ": number of station IDs does not match");
  }
  Data<double> lats;
  Data<double> lons;
  Data<double> julianDays;
  Data<double> depths;
  Data<std::string> stationTypes;
  Data<std::string> stationIDs;
  std::string juldReference;
  const size_t nObs;
  const size_t nLevels;
  const size_t nLocations;
};

/// \brief Interface to the NetCDF library to write feedback files
template <typename ncVarType = double>
class Writer {
 public:
  static std::string className()
    { return "nemo_feedback::feedback_io::Writer"; }
  Writer(
      eckit::PathName& filename,
      const MetaData & metaData,
      const NameData & name_data,
      const std::vector<bool> & isExtraVariable);

  /// \brief Write variable data
  void write_variable(
      const std::string & variable_name,
      const Data<ncVarType>& data);

  /// \brief Write surface variable data
  void write_variable_surf(
      const std::string & variable_name,
      const Data<ncVarType>& data);

  /// \brief Write profile variable data
  void write_variable_profile(
      const std::string & variable_name,
      const Data<ncVarType>& data);

  /// \brief Write surface QC data variable
  void write_variable_surf_qc(
      const std::string & variable_name,
      const Data<int32_t>& data);

  void write_variable_surf_qc(
      const std::string & variable_name,
      const Data<QC::Level>& data);

  /// \brief Write surface QC data variable with specified flag
  void write_variable_surf_qc(
      const std::string & variable_name,
      const Data<int32_t>& data,
      const size_t flag_index);

  /// \brief Write level QC data variable
  void write_variable_level_qc(
      const std::string & variable_name,
      const Data<int32_t>& data);

  void write_variable_level_qc(
      const std::string & variable_name,
      const Data<QC::Level>& data);

  /// \brief Write level QC data variable with specified flag
  void write_variable_level_qc(
      const std::string & variable_name,
      const Data<int32_t>& data,
      const size_t flag_index);

  static const netCDF::NcType typeToNcType() {
    return feedback_io::NetCDFTypeMap<ncVarType>::ncType;
  }

 private:
  Writer() : ncFile(), nobs_dim(), nlevels_dim(), metaData_(),
                         name_data_() {}

  /// \brief Define the coordinate variables in the NetCDF file
  void define_coord_variables(
      const size_t nObsVariables,
      const size_t nAdditionalEntries,
      const size_t nExtraVariables);

  /// \brief Write the coordinate variables to the NetCDF file
  void write_coord_variables();

  /// \brief Write the required feedback file metadata to the NetCDF file
  void write_metadata_variables(
      const std::vector<bool>& isExtraVariable);

  /// \brief Define the variables that impact entire observations in the
  ///        NetCDF file
  void define_whole_report_variables();

  /// \brief Write the variables that impact entire observations in the
  ///        NetCDF file
  void write_whole_report_variables();

  /// \brief Define a data variable in the NetCDF file
  void define_variable(
      const std::string & variable_name,
      const std::string & long_name,
      const std::string & unit_name,
      bool legacy_ops_qc_conventions = true);

  /// \brief Define an 'extra' data variable in the NetCDF file
  void define_extra_variable(
      const std::string & variable_name,
      const std::string & long_name,
      const std::string & unit_name,
      bool legacy_ops_qc_conventions = true);

  std::unique_ptr<netCDF::NcFile> ncFile;
  std::unique_ptr<netCDF::NcDim> nobs_dim;
  std::unique_ptr<netCDF::NcDim> nlevels_dim;
  std::unique_ptr<netCDF::NcDim> nqcf_dim;
  const MetaData metaData_;
  const NameData name_data_;
  static const std::map<std::string, size_t> coord_sizes;
};
}  // namespace feedback_io
}  // namespace nemo_feedback

//------------------------------------------------------------------------------
