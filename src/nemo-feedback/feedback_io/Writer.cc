/*
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "nemo-feedback/feedback_io/Writer.h"

#include <netcdf>
// Using Lynton Appel's netcdf-cxx4 from
// https://github.com/Unidata/netcdf-cxx4

#include <algorithm>
#include <sstream>
#include <iomanip>
#include <cstring>

#include "eckit/exception/Exceptions.h"

#include "oops/util/Logger.h"
#include "oops/util/Duration.h"

#include "nemo-feedback/feedback_io/Utils.h"

// names of netCDF dimensions, variables, attributes
#define N_QCF "N_QCF"
#define N_ENTRIES "N_ENTRIES"
#define N_EXTRA "N_EXTRA"
#define N_VARS "N_VARS"
#define N_OBS "N_OBS"
#define N_LEVELS "N_LEVELS"
#define STRINGNAM "STRINGNAM"
#define STRINGNAM_NUM 8
#define STRINGGRID "STRINGGRID"
#define STRINGWMO "STRINGWMO"
#define STRINGWMO_NUM 8
#define STRINGTYP "STRINGTYP"
#define STRINGTYP_NUM 4
#define STRINGJULD "STRINGJULD"
#define STRINGJULD_NUM 14

#define QC_CONVENTIONS "U.S. Integrated Ocean Observing System, 2017. Manual "\
  "for the Use of Real-Time Oceanographic Data Quality Control Flags, Version"\
  " 1.1"

namespace nemo_feedback {
namespace feedback_io {

template<class T>
const std::map<std::string, size_t> Writer<T>::coord_sizes{
  {N_QCF, 2}, {STRINGNAM, STRINGNAM_NUM}, {STRINGGRID, 1},
  {STRINGWMO, STRINGWMO_NUM}, {STRINGTYP, STRINGTYP_NUM},
  {STRINGJULD, STRINGJULD_NUM}};

typedef char fixed_length_name_type[STRINGNAM_NUM+1];
typedef char fixed_length_type_type[STRINGTYP_NUM+1];

template<class C>
Writer<C>::Writer(
      eckit::PathName& filename,
      const MetaData & metaData,
      const NameData & name_data,
      const std::vector<bool> & isExtraVariable)
    : ncFile(nullptr), metaData_(metaData),
      name_data_(name_data) {
  oops::Log::trace() << "nemo_feedback::NemoFieldWriter::NemoFieldWriter"
                     << " constructing for: " << filename.fullName().asString()
                     << std::endl;

  ncFile = std::make_unique<netCDF::NcFile>(filename.fullName().asString(),
      netCDF::NcFile::replace);
  if (ncFile->isNull()) {
    std::ostringstream err_stream;
    err_stream << "nemo_feedback::NemoFieldWriter::NemoFieldWriter cannot open "
               << filename << std::endl;
    throw eckit::BadValue(err_stream.str(), Here());
  }

  const size_t nExtraVariables = std::count(isExtraVariable.begin(),
      isExtraVariable.end(), true);
  const size_t nObsVariables = name_data_.variable_names.size()
    - nExtraVariables;
  const size_t nAdditionalEntries = name_data_.additional_names.size();

  name_data_.validate();

  define_coord_variables(
      nObsVariables,
      nAdditionalEntries,
      nExtraVariables);
  write_metadata_variables(isExtraVariable);
  write_coord_variables();
  define_whole_report_variables();
  if (metaData_.nObs != 0) {
    write_whole_report_variables();
  }

  for (size_t i = 0; i < name_data_.variable_names.size(); ++i) {
    if (isExtraVariable[i]) {
      define_extra_variable(name_data_.variable_names[i],
                            name_data_.long_names[i],
                            name_data_.unit_names[i],
                            name_data_.legacy_ops_qc_conventions[i]);
    } else {
      define_variable(
          name_data_.variable_names[i],
          name_data_.long_names[i],
          name_data_.unit_names[i],
          name_data_.legacy_ops_qc_conventions[i]);
    }
  }
}

template <class C>
void Writer<C>::define_coord_variables(
    const size_t nObsVariables,
    const size_t nAdditionalEntries,
    const size_t nExtraVariables) {
  netCDF::NcDim tmp_Dim;
  for (const auto& kv : coord_sizes) {
    if (kv.first == N_QCF) {
      nqcf_dim = std::make_unique<netCDF::NcDim>(
          ncFile->addDim(kv.first, kv.second));
    } else {
      tmp_Dim = ncFile->addDim(kv.first, kv.second);
    }
  }

  nobs_dim = std::make_unique<netCDF::NcDim>(
      ncFile->addDim(N_OBS, metaData_.nObs));
  nlevels_dim = std::make_unique<netCDF::NcDim>(ncFile->addDim(N_LEVELS,
        metaData_.nLevels));
  tmp_Dim = ncFile->addDim(N_VARS, nObsVariables);
  tmp_Dim = ncFile->addDim(N_ENTRIES, nAdditionalEntries);
  if (nExtraVariables > 0) {
    tmp_Dim = ncFile->addDim(N_EXTRA, nExtraVariables);
  }
}

template <class C>
void Writer<C>::write_metadata_variables(
    const std::vector<bool>& isExtraVariable) {
  {
    ncFile->putAtt("title", "NEMO observation operator output");
    ncFile->putAtt("Convention", "NEMO unified observation operator output");
    netCDF::NcVar nc_juld_var = ncFile->addVar("JULD_REFERENCE", netCDF::ncChar,
        ncFile->getDim(STRINGJULD));
    nc_juld_var.putAtt("long_name", "Date of reference for julian days");
    nc_juld_var.putAtt("Conventions", "YYYYMMDDHHMMSS");

    char buffer[15];
    metaData_.juldReference.copy(buffer, 14);
    buffer[14] = '\0';
    nc_juld_var.putVar({0}, {STRINGJULD_NUM}, buffer);
  }

  {
    size_t n_vars = name_data_.variable_names.size();
    int nExtraVariables = std::count(isExtraVariable.begin(),
        isExtraVariable.end(), true);
    n_vars -= nExtraVariables;

    std::vector<netCDF::NcDim>
        dims{ncFile->getDim(N_VARS), ncFile->getDim(STRINGNAM)};
    netCDF::NcVar nc_var_list_var = ncFile->addVar("VARIABLES",
        netCDF::ncChar, dims);
    netCDF::NcVar nc_var_list_extra;
    if (nExtraVariables > 0) {
      std::vector<netCDF::NcDim>
          dimsextra{ncFile->getDim(N_EXTRA), ncFile->getDim(STRINGNAM)};
      nc_var_list_extra = ncFile->addVar("EXTRA",
                                         netCDF::ncChar,
                                         dimsextra);
    }
    nc_var_list_var.putAtt("long_name", "List of variables in feedback files");
    fixed_length_name_type data;
    size_t ivar = 0;
    size_t iextra = 0;
    // trim and pad string to fit in nc char array
    for (size_t i=0; i < (n_vars + nExtraVariables); ++i) {
      for (size_t j=0; j < STRINGNAM_NUM; ++j) {
        if (j < name_data_.variable_names[i].length()) {
          data[j] = static_cast<char>(name_data_.variable_names.at(i).at(j));
        } else {data[j] = ' ';}
      }
      if (isExtraVariable[i]) {
        nc_var_list_extra.putVar({iextra++, 0}, {1, STRINGNAM_NUM}, data);
      } else {
        nc_var_list_var.putVar({ivar++, 0}, {1, STRINGNAM_NUM}, data);
      }
    }
  }

  {
    size_t nAdditionalEntries = name_data_.additional_names.size();
    std::vector<netCDF::NcDim> dims{ncFile->getDim(N_ENTRIES),
        ncFile->getDim(STRINGNAM)};
    netCDF::NcVar nc_entries_var = ncFile->addVar("ENTRIES",
        netCDF::ncChar, dims);
    nc_entries_var.putAtt("long_name",
        "List of additional entries for each variable in feedback files");
    fixed_length_name_type data;;
    // trim and pad string to fit in nc char array
    for (size_t i=0; i < nAdditionalEntries; ++i) {
      for (size_t j=0; j < STRINGNAM_NUM; ++j) {
        if (j < name_data_.additional_names[i].length()) {
          data[j] = static_cast<char>(name_data_.additional_names.at(i).at(j));
        } else {data[j] = ' ';}
      }
      nc_entries_var.putVar({i, 0}, {1, STRINGNAM_NUM}, data);
    }
  }
}

/// define metadata variables for each report, applying to all observations in
/// that report. e.g: STATION_IDENTIFIER, STATION_TYPE, DEPTH_QC,
/// DEPTH_QC_FLAGS, OBSERVATION_QC, OBSERVATION_QC_FLAGS, JULD_QC,
/// JULD_QC_FLAGS, POSITION_QC, POSITION_QC_FLAGS
template <class C>
void Writer<C>::define_whole_report_variables() {
  std::string flag_conventions = "JEDI UFO QC flag conventions";
  {
    const std::vector<netCDF::NcDim> dims{*nobs_dim, ncFile->getDim(STRINGWMO)};
    netCDF::NcVar tmp_var = ncFile->addVar("STATION_IDENTIFIER", netCDF::ncChar,
        dims);
    tmp_var.putAtt("long_name", "Station identifier");
  }

  {
    const std::vector<netCDF::NcDim> dims{*nobs_dim, ncFile->getDim(STRINGTYP)};
    netCDF::NcVar tmp_var = ncFile->addVar("STATION_TYPE", netCDF::ncChar,
        dims);
    tmp_var.putAtt("long_name", "Code instrument type");
  }

  {
    const std::vector<netCDF::NcDim> dims{*nobs_dim, *nlevels_dim};
    netCDF::NcVar tmp_var = ncFile->addVar("DEPTH_QC", netCDF::ncInt, dims);
    tmp_var.putAtt("long_name", "Quality on depth");
    tmp_var.putAtt("Conventions", QC_CONVENTIONS);
  }

  {
    const std::vector<netCDF::NcDim> dims{*nobs_dim, *nlevels_dim,
        ncFile->getDim(N_QCF)};
    netCDF::NcVar tmp_var = ncFile->addVar("DEPTH_QC_FLAGS", netCDF::ncInt,
        dims);
    tmp_var.putAtt("long_name", "Quality on depth");
    tmp_var.putAtt("Conventions", flag_conventions);
  }

  {
    netCDF::NcVar tmp_var = ncFile->addVar("OBSERVATION_QC", netCDF::ncInt,
        *nobs_dim);
    tmp_var.putAtt("long_name", "Quality on observation");
    tmp_var.putAtt("Conventions", QC_CONVENTIONS);
  }

  {
    const std::vector<netCDF::NcDim> dims{*nobs_dim, ncFile->getDim(N_QCF)};
    netCDF::NcVar tmp_var = ncFile->addVar("OBSERVATION_QC_FLAGS",
        netCDF::ncInt, dims);
    tmp_var.putAtt("long_name", "Quality on observation");
    tmp_var.putAtt("Conventions", flag_conventions);
  }

  {
    netCDF::NcVar tmp_var = ncFile->addVar("POSITION_QC", netCDF::ncInt,
        *nobs_dim);
    tmp_var.putAtt("long_name", "Quality on position");
    tmp_var.putAtt("Conventions", QC_CONVENTIONS);
  }

  {
    const std::vector<netCDF::NcDim> dims{*nobs_dim, ncFile->getDim(N_QCF)};
    netCDF::NcVar tmp_var = ncFile->addVar("POSITION_QC_FLAGS", netCDF::ncInt,
        dims);
    tmp_var.putAtt("long_name", "Quality on position");
    tmp_var.putAtt("Conventions", flag_conventions);
  }

  {
    netCDF::NcVar tmp_var = ncFile->addVar("JULD_QC", netCDF::ncInt, *nobs_dim);
    tmp_var.putAtt("long_name", "Quality on date and time");
    tmp_var.putAtt("Conventions", QC_CONVENTIONS);
  }

  {
    const std::vector<netCDF::NcDim> dims{*nobs_dim, ncFile->getDim(N_QCF)};
    netCDF::NcVar tmp_var = ncFile->addVar("JULD_QC_FLAGS", netCDF::ncInt,
        dims);
    tmp_var.putAtt("long_name", "Quality on date and time");
    tmp_var.putAtt("Conventions", flag_conventions);
  }

  {
    netCDF::NcVar tmp_var = ncFile->addVar("ORIGINAL_FILE_INDEX", netCDF::ncInt,
        *nobs_dim);
    tmp_var.putAtt("long_name", "Index in original data file");
  }
}

/// Define a main variable with its auxillary variables in the feedback file
/// e.g: SST_ +
/// "OBS", "Hx", "QC", "QC_FLAGS", "LEVEL_QC",
/// "LEVEL_QC_FLAGS", "IOBSI", "IOBSJ", "IOBSK"
/// and additional named variables passed in to the method. All variables are
/// assumed to be on the T grid.
template <class C>
void Writer<C>::define_variable(
    const std::string & variable_name,
    const std::string & longName,
    const std::string & units,
    bool legacy_ops_qc_conventions) {

  std::string flag_conventions = "JEDI UFO QC flag conventions";
  if (legacy_ops_qc_conventions)
    flag_conventions = "OPS flag conventions";

  const std::vector<netCDF::NcDim> dims{*nobs_dim, *nlevels_dim};
  netCDF::NcVar obs_var;
  obs_var = ncFile->addVar(variable_name + "_OBS",
            typeToNcType(), dims);
  obs_var.putAtt("long_name", longName);
  obs_var.putAtt("units", units);
  obs_var.setFill(true, feedback_io::typeToFill::value<C>());

  for (const auto &name : name_data_.additional_names) {
    netCDF::NcVar add_var;
    add_var = ncFile->addVar(variable_name + "_" + name,
            typeToNcType(), dims);
    add_var.putAtt("long_name", longName + " " + name);
    add_var.putAtt("units", units);
    add_var.setFill(true, feedback_io::typeToFill::value<C>());
  }

  {
    const std::vector<netCDF::NcDim> qcf_dims{*nobs_dim, *nqcf_dim};
    netCDF::NcVar qc_flags_var = ncFile->addVar(variable_name + "_QC_FLAGS",
        netCDF::ncInt, qcf_dims);
    qc_flags_var.setFill(true, feedback_io::typeToFill::value<int32_t>());
    qc_flags_var.putAtt("long_name", std::string("quality flags on ")
                        + longName);
    qc_flags_var.putAtt("Conventions", flag_conventions);
  }

  {
    const std::vector<netCDF::NcDim> lvl_qcf_dims{*nobs_dim, *nlevels_dim,
        *nqcf_dim};
    netCDF::NcVar level_qc_flags_var = ncFile->addVar(variable_name
        + "_LEVEL_QC_FLAGS", netCDF::ncInt, lvl_qcf_dims);
    level_qc_flags_var.setFill(true, feedback_io::typeToFill::value<int32_t>());
    level_qc_flags_var.putAtt("long_name",
        std::string("quality flags for each level on ") + longName);
    level_qc_flags_var.putAtt("Conventions", flag_conventions);
  }

  netCDF::NcVar qc_var = ncFile->addVar(variable_name + "_QC", netCDF::ncInt,
      *nobs_dim);
  qc_var.putAtt("long_name", std::string("quality on ") + longName);
  qc_var.putAtt("Conventions", QC_CONVENTIONS);

  netCDF::NcVar level_qc_var = ncFile->addVar(variable_name + "_LEVEL_QC",
      netCDF::ncInt, dims);
  level_qc_var.setFill(true, feedback_io::typeToFill::value<int32_t>());
  level_qc_var.putAtt("long_name", std::string("quality for each level on ")
                      + longName);
  level_qc_var.putAtt("Conventions", QC_CONVENTIONS);

  {
    netCDF::NcVar tmp_var;
    tmp_var = ncFile->addVar(variable_name + "_IOBSI", netCDF::ncInt,
        *nobs_dim);
    tmp_var.putAtt("long_name", "ORCA grid search I coordinate");
    tmp_var = ncFile->addVar(variable_name + "_IOBSJ", netCDF::ncInt,
        *nobs_dim);
    tmp_var.putAtt("long_name", "ORCA grid search J coordinate");
    tmp_var = ncFile->addVar(variable_name + "_IOBSK", netCDF::ncInt, dims);
    tmp_var.putAtt("long_name", "ORCA grid search K coordinate");
    tmp_var = ncFile->addVar(variable_name + "_GRID", netCDF::ncChar,
        ncFile->getDim(STRINGGRID));
    tmp_var.putAtt("long_name", "ORCA grid search grid (T,U,V)");
    char data[2] {"T"};
    tmp_var.putVar({0}, {1}, data);
  }
}

/// Define an extra variable.
template <class C>
void Writer<C>::define_extra_variable(
    const std::string & variable_name,
    const std::string & longName,
    const std::string & units,
    bool legacy_ops_qc_conventions) {

  std::string flag_conventions = "JEDI UFO QC flag conventions";
  if (legacy_ops_qc_conventions)
    flag_conventions = "OPS flag conventions";

  const std::vector<netCDF::NcDim> dims{*nobs_dim, *nlevels_dim};

  netCDF::NcVar var = ncFile->addVar(variable_name,
                                     netCDF::ncDouble, dims);
  var.putAtt("long_name", longName);
  var.putAtt("units", units);
}

template <class C>
void Writer<C>::write_coord_variables() {
  netCDF::NcVar lat_var = ncFile->addVar("LATITUDE", netCDF::ncDouble,
      *nobs_dim);
  lat_var.putAtt("units", "degrees_north");
  lat_var.putAtt("long_name", "latitude");

  netCDF::NcVar lon_var = ncFile->addVar("LONGITUDE", netCDF::ncDouble,
      *nobs_dim);
  lon_var.putAtt("units", "degrees_east");
  lon_var.putAtt("long_name", "longitude");

  const std::vector<netCDF::NcDim> dims{*nobs_dim, *nlevels_dim};
  netCDF::NcVar depth_var = ncFile->addVar("DEPTH", netCDF::ncDouble,
      dims);
  depth_var.putAtt("units", "metre");
  depth_var.putAtt("long_name", "Depth");

  netCDF::NcVar juld_var = ncFile->addVar("JULD", netCDF::ncDouble, *nobs_dim);
  juld_var.putAtt("units", "days since JULD_REFERENCE");
  juld_var.putAtt("long_name", "Julian day");

  if (metaData_.nObs == 0) return;

  std::vector<double> buffer = metaData_.lats.raw_surface();
  lat_var.putVar(buffer.data());

  buffer = metaData_.lons.raw_surface();
  lon_var.putVar(buffer.data());

  for (size_t iProfile = 0; iProfile < metaData_.depths.n_obs(); ++iProfile) {
    auto profileBuffer = metaData_.depths.raw_profile(iProfile);
    depth_var.putVar({iProfile, 0},
                     {1, profileBuffer.size()},
                     profileBuffer.data());
  }

  buffer = metaData_.julianDays.raw_surface();
  juld_var.putVar(buffer.data());
}

template <class C>
void Writer<C>::write_whole_report_variables() {
  // Write station type.
  {
    size_t nchars = (ncFile->getDim(STRINGTYP)).getSize();
    auto station_type_var = ncFile->getVar("STATION_TYPE");
    // +1 is for the null-terminator of a cstring
    std::vector<std::string> stationTypesSurface = metaData_.stationTypes
      .raw_surface();
    ASSERT_MSG(stationTypesSurface.size() == metaData_.nObs,
        "stationTypesSurface.size() != metaData_.nObs");
    std::string buffer;
    buffer.reserve(metaData_.nObs*nchars);
    for (std::string stationType : stationTypesSurface) {
      stationType.resize(STRINGTYP_NUM, ' ');
      buffer += stationType.substr(0, STRINGTYP_NUM);
    }
    station_type_var.putVar({0, 0},
                            {metaData_.nObs, nchars},
                            buffer.c_str());
  }

  // Write station IDs.
  {
    size_t nchars = (ncFile->getDim(STRINGWMO)).getSize();
    auto station_id_var = ncFile->getVar("STATION_IDENTIFIER");
    // +1 is for the null-terminator of a cstring
    std::vector<std::string> stationIDsSurface = metaData_.stationIDs
      .raw_surface();
    ASSERT_MSG(stationIDsSurface.size() == metaData_.nObs,
        "stationTypesSurface.size() != metaData_.nObs");
    std::string buffer;
    buffer.reserve(metaData_.nObs*nchars);
    for (std::string stationID : stationIDsSurface) {
      stationID.resize(STRINGWMO_NUM, ' ');
      buffer += stationID.substr(0, STRINGWMO_NUM);
    }
    station_id_var.putVar({0, 0},
                          {metaData_.nObs, nchars},
                          buffer.c_str());
  }
}

template <class C>
void Writer<C>::write_variable(
    const std::string & variable_name,
    const Data<C>& data) {
  if (metaData_.nLevels <= 1) {
    write_variable_surf(variable_name, data);
  } else {
    write_variable_profile(variable_name, data);
  }
}

template <class C>
void Writer<C>::write_variable_surf(
    const std::string & variable_name,
    const Data<C>& data) {
  oops::Log::trace() << Writer::className() << "::write_variable_surf: writing "
                     << variable_name << std::endl;
  auto surf_var = ncFile->getVar(variable_name);
  if (data.n_locations() == 0) return;
  std::vector<C> buffer = data.raw_surface();
  surf_var.putVar(buffer.data());
}

template <class C>
void Writer<C>::write_variable_surf_qc(
    const std::string & variable_name,
    const Data<int32_t>& data) {
  oops::Log::trace() << Writer::className()
                     << "::write_variable_surf_qc: writing "
                     << variable_name << std::endl;
  auto surf_var = ncFile->getVar(variable_name);
  if (data.n_locations() == 0) return;
  std::vector<int32_t> buffer = data.raw_surface();
  surf_var.putVar(buffer.data());
}

template <class C>
void Writer<C>::write_variable_surf_qc(
    const std::string & variable_name,
    const Data<QC::Level>& data) {
  oops::Log::trace() << Writer::className()
                     << "::write_variable_surf_qc: writing "
                     << variable_name << std::endl;
  auto surf_var = ncFile->getVar(variable_name);
  if (data.n_locations() == 0) return;
  std::vector<QC::Level> buffer = data.raw_surface();
  surf_var.putVar(buffer.data());
}

template <class C>
void Writer<C>::write_variable_surf_qc(
    const std::string& variable_name,
    const Data<int32_t>& data,
    const size_t flag_index) {
  oops::Log::trace() << Writer::className()
                     << "::write_variable_surf_qc: writing flag_index: "
                     << flag_index << " of " << variable_name << std::endl;
  auto surf_var = ncFile->getVar(variable_name);
  if (data.n_locations() == 0) return;
  std::vector<int32_t> buffer = data.raw_surface();
  surf_var.putVar({0, flag_index}, {metaData_.nObs, 1}, buffer.data());
}

template <class C>
void Writer<C>::write_variable_profile(
    const std::string& variable_name,
    const Data<C>& data) {
  oops::Log::trace() << Writer::className()
                     << "::write_variable_profile: writing "
                     << variable_name
                     << " nLevels: " << metaData_.nLevels << std::endl;
  auto var = ncFile->getVar(variable_name);
  if (var.isNull()) {
      std::ostringstream err_stream;
      err_stream << Writer::className() << "::write_variable_profile "
                 << "ncVar '" << variable_name << "' is not present in "
                 <<"NetCDF file";
      throw eckit::BadValue(err_stream.str(), Here());
  }

  if (data.n_locations() == 0) return;

  if (metaData_.nLocations != data.n_locations()) {
      std::ostringstream err_stream;
      err_stream << Writer::className() << "::write_variable_profile "
                 << "index range out of bounds '" << variable_name << "' "
                 << metaData_.nLocations << " != " << data.n_locations();
      throw eckit::BadValue(err_stream.str(), Here());
  }
  for (size_t iProfile = 0; iProfile < data.n_obs(); ++iProfile) {
    auto profileBuffer = data.raw_profile(iProfile);
    var.putVar({iProfile, 0},
               {1, profileBuffer.size()},
               profileBuffer.data());
  }
}

template <class C>
void Writer<C>::write_variable_level_qc(
    const std::string & variable_name,
    const Data<int32_t>& data) {
  oops::Log::trace() << Writer::className()
                     << "::write_variable_level_qc: writing "
                     << variable_name << std::endl;
  auto var = ncFile->getVar(variable_name);
  if (var.isNull()) {
      std::ostringstream err_stream;
      err_stream << Writer::className() << "::write_variable_level_qc"
                 << " ncVar '" << variable_name << "' is not present in "
                 <<"NetCDF file";
      throw eckit::BadValue(err_stream.str(), Here());
    }

  if (data.n_locations() == 0) return;

  if (metaData_.nLocations != data.n_locations()) {
      std::ostringstream err_stream;
      err_stream << Writer::className() << "::write_variable_level_qc "
                 << "index range out of bounds '" << variable_name << "' "
                 << metaData_.nLocations << " != " << data.n_locations();
      throw eckit::BadValue(err_stream.str(), Here());
  }
  for (size_t iProfile = 0; iProfile < data.n_obs(); ++iProfile) {
    std::vector<int32_t> profileBuffer = data.raw_profile(iProfile);
    var.putVar({iProfile, 0},
               {1, profileBuffer.size()},
               profileBuffer.data());
  }
}

template <class C>
void Writer<C>::write_variable_level_qc(
    const std::string & variable_name,
    const Data<QC::Level>& data) {
  oops::Log::trace() << Writer::className()
                     << "::write_variable_level_qc: writing "
                     << variable_name << std::endl;
  auto var = ncFile->getVar(variable_name);
  if (var.isNull()) {
      std::ostringstream err_stream;
      err_stream << Writer::className() << "::write_variable_level_qc"
                 << " ncVar '" << variable_name << "' is not present in "
                 <<"NetCDF file";
      throw eckit::BadValue(err_stream.str(), Here());
    }

  if (data.n_locations() == 0) return;

  if (metaData_.nLocations != data.n_locations()) {
      std::ostringstream err_stream;
      err_stream << Writer::className() << "::write_variable_level_qc "
                 << "index range out of bounds '" << variable_name << "' "
                 << metaData_.nLocations << " != " << data.n_locations();
      throw eckit::BadValue(err_stream.str(), Here());
  }
  for (size_t iProfile = 0; iProfile < data.n_obs(); ++iProfile) {
    std::vector<QC::Level> buffer = data.raw_profile(iProfile);
    std::vector<int32_t> profileBuffer;
    std::transform(buffer.begin(), buffer.end(), profileBuffer.begin(),
        [](QC::Level value) {return static_cast<int32_t>(value);});
    var.putVar({iProfile, 0},
               {1, profileBuffer.size()},
               profileBuffer.data());
  }
}

template <class C>
void Writer<C>::write_variable_level_qc(
    const std::string & variable_name,
    const Data<int32_t>& data,
    const size_t flag_index) {
  oops::Log::trace() << Writer::className()
                     << "::write_variable_level_qc: writing flag_index: "
                     << flag_index << " of " << variable_name
                     << std::endl;
  auto var = ncFile->getVar(variable_name);
  if (var.isNull()) {
      std::ostringstream err_stream;
      err_stream << Writer::className() << "::write_variable_level_qc"
                 << " ncVar '" << variable_name << "' is not present in "
                 <<"NetCDF file";
      throw eckit::BadValue(err_stream.str(), Here());
    }

  if (data.n_locations() == 0) return;

  if (metaData_.nLocations != data.n_locations()) {
      std::ostringstream err_stream;
      err_stream << Writer::className() << "::write_variable_level_qc "
                 << "index range out of bounds '" << variable_name << "' "
                 << metaData_.nLocations << " != " << data.n_locations();
      throw eckit::BadValue(err_stream.str(), Here());
  }
  for (size_t iProfile = 0; iProfile < data.n_obs(); ++iProfile) {
    auto profileBuffer = data.raw_profile(iProfile);
    var.putVar({iProfile, 0, flag_index},
               {1, profileBuffer.size(), 1},
               profileBuffer.data());
  }
}

template class Writer<double>;
template class Writer<float>;

}  // namespace feedback_io
}  // namespace nemo_feedback
