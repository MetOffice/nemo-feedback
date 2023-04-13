/*
 * (C) British Crown Copyright 2021 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "nemo-feedback/NemoFeedbackWriter.h"

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

// default
template <typename T> const netCDF::NcType NetCDFTypeMap<T>::ncType =
  netCDF::ncDouble;
template <> const netCDF::NcType NetCDFTypeMap<float>::ncType = netCDF::ncFloat;

template <> const double typeToFill::value<double>() {
  return NemoFeedbackWriter<double>::double_fillvalue;
}
template <> const float  typeToFill::value<float>() {
  return NemoFeedbackWriter<double>::float_fillvalue;
}
template <> const int32_t  typeToFill::value<int32_t>() {
  return NemoFeedbackWriter<double>::int32_fillvalue;
}

template<class T>
const std::map<std::string, size_t> NemoFeedbackWriter<T>::coord_sizes{
  {N_QCF, 2}, {STRINGNAM, STRINGNAM_NUM}, {STRINGGRID, 1},
  {STRINGWMO, STRINGWMO_NUM}, {STRINGTYP, STRINGTYP_NUM},
  {STRINGJULD, STRINGJULD_NUM}};

typedef char fixed_length_name_type[STRINGNAM_NUM+1];
typedef char fixed_length_type_type[STRINGTYP_NUM+1];

void NameData::validate() const {
  const size_t n_vars = variable_names.size();
  if ((legacy_ops_qc_conventions.size() != n_vars) ||
      (long_names.size() != n_vars) ||
      (unit_names.size() != n_vars)) {
    std::ostringstream err_stream;
    err_stream << "nemo_feedback::NameData size mismatch " << std::endl;
    throw eckit::BadValue(err_stream.str(), Here());
  }
}

template<class C>
NemoFeedbackWriter<C>::NemoFeedbackWriter(
    eckit::PathName& filename,
    const CoordData & coords,
    const NameData & name_data,
    const std::vector<bool> & extra_vars,
    const std::vector<std::string> & station_types,
    const std::vector<std::string> & station_ids)
    : ncFile(nullptr), coords_(coords),
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

  const size_t n_extra = std::count(extra_vars.begin(), extra_vars.end(), true);
  const size_t n_obs_vars = name_data_.variable_names.size() - n_extra;
  const size_t max_n_add_entries = name_data_.additional_names.size();

  name_data_.validate();

  define_coord_variables(
      n_obs_vars,
      max_n_add_entries,
      n_extra);
  write_metadata_variables(extra_vars);
  write_coord_variables();
  define_whole_report_variables();
  if (coords_.n_obs != 0) {
    write_whole_report_variables(station_types, station_ids);
  }

  for (int i=0; i < name_data_.variable_names.size(); ++i) {
    if (extra_vars[i]) {
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
void NemoFeedbackWriter<C>::define_coord_variables(
    const size_t n_obs_vars,
    const size_t n_add_entries,
    const size_t n_extra) {
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
      ncFile->addDim(N_OBS, coords_.n_obs));
  nlevels_dim = std::make_unique<netCDF::NcDim>(ncFile->addDim(N_LEVELS,
        coords_.n_levels));
  tmp_Dim = ncFile->addDim(N_VARS, n_obs_vars);
  tmp_Dim = ncFile->addDim(N_ENTRIES, n_add_entries);
  if (n_extra > 0) {
    tmp_Dim = ncFile->addDim(N_EXTRA, n_extra);
  }
}

template <class C>
void NemoFeedbackWriter<C>::write_metadata_variables(
    const std::vector<bool>& extra_vars) {
  {
    ncFile->putAtt("title", "NEMO observation operator output");
    ncFile->putAtt("Convention", "NEMO unified observation operator output");
    netCDF::NcVar nc_juld_var = ncFile->addVar("JULD_REFERENCE", netCDF::ncChar,
        ncFile->getDim(STRINGJULD));
    nc_juld_var.putAtt("long_name", "Date of reference for julian days");
    nc_juld_var.putAtt("Conventions", "YYYYMMDDHHMMSS");

    int year, month, day, hour, minute, second;
    coords_.juld_reference.toYYYYMMDDhhmmss(year, month, day, hour, minute,
                                           second);
    std::ostringstream ref_stream;
    ref_stream << std::setfill('0');
    ref_stream << std::setw(4) << year;
    ref_stream << std::setw(2) << month;
    ref_stream << std::setw(2) << day;
    ref_stream << std::setw(2) << hour;
    ref_stream << std::setw(2) << minute;
    ref_stream << std::setw(2) << second;

    char buffer[15];
    ref_stream.str().copy(buffer, 14);
    buffer[14] = '\0';
    nc_juld_var.putVar({0}, {STRINGJULD_NUM}, buffer);
  }

  {
    size_t n_vars = name_data_.variable_names.size();
    int n_extra = std::count(extra_vars.begin(), extra_vars.end(), true);
    n_vars -= n_extra;

    std::vector<netCDF::NcDim>
        dims{ncFile->getDim(N_VARS), ncFile->getDim(STRINGNAM)};
    netCDF::NcVar nc_var_list_var = ncFile->addVar("VARIABLES",
        netCDF::ncChar, dims);
    netCDF::NcVar nc_var_list_extra;
    if (n_extra > 0) {
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
    for (size_t i=0; i < (n_vars + n_extra); ++i) {
      for (size_t j=0; j < STRINGNAM_NUM; ++j) {
        if (j < name_data_.variable_names[i].length()) {
          data[j] = static_cast<char>(name_data_.variable_names.at(i).at(j));
        } else {data[j] = ' ';}
      }
      if (extra_vars[i]) {
        nc_var_list_extra.putVar({iextra++, 0}, {1, STRINGNAM_NUM}, data);
      } else {
        nc_var_list_var.putVar({ivar++, 0}, {1, STRINGNAM_NUM}, data);
      }
    }
  }

  {
    size_t max_n_add_entries = name_data_.additional_names.size();
    std::vector<netCDF::NcDim> dims{ncFile->getDim(N_ENTRIES),
        ncFile->getDim(STRINGNAM)};
    netCDF::NcVar nc_entries_var = ncFile->addVar("ENTRIES",
        netCDF::ncChar, dims);
    nc_entries_var.putAtt("long_name",
        "List of additional entries for each variable in feedback files");
    fixed_length_name_type data;;
    // trim and pad string to fit in nc char array
    for (size_t i=0; i < max_n_add_entries; ++i) {
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
void NemoFeedbackWriter<C>::define_whole_report_variables() {

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
void NemoFeedbackWriter<C>::define_variable(
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
  obs_var.setFill(true, typeToFill::value<C>());

  for (const auto &name : name_data_.additional_names) {
    netCDF::NcVar add_var;
    add_var = ncFile->addVar(variable_name + "_" + name,
            typeToNcType(), dims);
    add_var.putAtt("long_name", longName + " " + name);
    add_var.putAtt("units", units);
    add_var.setFill(true, typeToFill::value<C>());
  }

  {
    const std::vector<netCDF::NcDim> qcf_dims{*nobs_dim, *nqcf_dim};
    netCDF::NcVar qc_flags_var = ncFile->addVar(variable_name + "_QC_FLAGS",
        netCDF::ncInt, qcf_dims);
    qc_flags_var.setFill(true, int32_fillvalue);
    qc_flags_var.putAtt("long_name", std::string("quality flags on ")
                        + longName);
    qc_flags_var.putAtt("Conventions", flag_conventions);
  }

  {
    const std::vector<netCDF::NcDim> lvl_qcf_dims{*nobs_dim, *nlevels_dim,
        *nqcf_dim};
    netCDF::NcVar level_qc_flags_var = ncFile->addVar(variable_name
        + "_LEVEL_QC_FLAGS", netCDF::ncInt, lvl_qcf_dims);
    level_qc_flags_var.setFill(true, int32_fillvalue);
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
  level_qc_var.setFill(true, int32_fillvalue);
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
void NemoFeedbackWriter<C>::define_extra_variable(
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
void NemoFeedbackWriter<C>::write_coord_variables() {
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

  if (coords_.n_obs == 0) return;

  lat_var.putVar(coords_.lats.data());

  lon_var.putVar(coords_.lons.data());

  if (coords_.n_obs == coords_.n_locs) {
    depth_var.putVar(coords_.depths.data());
  } else {
    for (size_t n = 0; n < coords_.n_obs; ++n) {
      depth_var.putVar({n, 0},
                       {1, coords_.record_counts[n]},
                       coords_.depths.data() + coords_.record_starts[n]);
    }
  }

  juld_var.putVar(coords_.julian_days.data());
}

template <class C>
void NemoFeedbackWriter<C>::write_whole_report_variables(
    const std::vector<std::string> & station_types,
    const std::vector<std::string> & station_ids) {

  if (station_ids.size() != coords_.n_obs ||
      station_types.size() != coords_.n_obs )
    throw eckit::BadValue(std::string("NemoFeedbackWriter::") +
        "write_whole_report_variables: station_ids or station_types " +
        "not defined for some observations", Here());

  // Write station type.
  {
    size_t nchars = (ncFile->getDim(STRINGTYP)).getSize();
    auto station_type_var = ncFile->getVar("STATION_TYPE");
    int j = 0;
    // +1 is for the null-terminator of a cstring
    char* buffer = new char[coords_.n_obs*nchars+1];
    for (int i = 0; i < coords_.n_obs; ++i) {
      if (coords_.n_obs == coords_.n_locs) {
        strcpy(buffer + nchars*j++, station_types[i].c_str());
      } else {
        strcpy(buffer + nchars*j++, station_types[i].c_str());
      }
    }
    station_type_var.putVar({0, 0},
                            {coords_.n_obs, nchars},
                            buffer);
    delete[] buffer;
  }

  // Write station IDs.
  {
    size_t nchars = (ncFile->getDim(STRINGWMO)).getSize();
    auto station_id_var = ncFile->getVar("STATION_IDENTIFIER");
    int j = 0;
    // +1 is for the null-terminator of a cstring
    char* buffer = new char[coords_.n_obs*nchars+1];
    for (int i = 0; i < coords_.n_obs; ++i) {
      strcpy(buffer + nchars*j++, station_ids[i].c_str());
    }
    station_id_var.putVar({0, 0},
                          {coords_.n_obs, nchars},
                          buffer);
    delete[] buffer;
  }
}

template <class C>
void NemoFeedbackWriter<C>::write_variable_surf(
    const std::string & variable_name,
    const std::vector<C>& data) {
  oops::Log::trace() << "NemoFeedbackWriter::write_variable_surf: writing "
                     << variable_name << std::endl;
  auto surf_var = ncFile->getVar(variable_name);
  if (data.size() == 0) return;
  surf_var.putVar(data.data());
}

template <class C>
void NemoFeedbackWriter<C>::write_variable_surf_qc(
    const std::string & variable_name,
    const std::vector<int32_t>& data) {
  oops::Log::trace() << "NemoFeedbackWriter::write_variable_surf_qc: writing "
                     << variable_name << std::endl;
  auto surf_var = ncFile->getVar(variable_name);
  if (data.size() == 0) return;
  surf_var.putVar(data.data());
}

template <class C>
void NemoFeedbackWriter<C>::write_variable_surf_qc(
    const std::string& variable_name,
    const std::vector<int32_t>& data,
    const size_t flag_index) {
  oops::Log::trace() << "NemoFeedbackWriter::write_variable_surf_qc: writing "
                     << "flag_index: " << flag_index << " of " << variable_name
                     << std::endl;
  auto surf_var = ncFile->getVar(variable_name);
  if (data.size() == 0) return;
  surf_var.putVar({0, flag_index}, {coords_.n_obs, 1}, data.data());
}

template <class C>
void NemoFeedbackWriter<C>::write_variable_profile(
    const std::string& variable_name,
    const std::vector<C>& data) {
  oops::Log::trace() << "NemoFeedbackWriter::write_variable_profile: writing "
                     << variable_name << std::endl;
  auto var = ncFile->getVar(variable_name);
  if (var.isNull()) {
      std::ostringstream err_stream;
      err_stream << "nemo_feedback::NemoFeedbackWriter::write_variable_profile "
                 << "ncVar '" << variable_name << "' is not present in "
                 <<"NetCDF file";
      throw eckit::BadValue(err_stream.str(), Here());
  }

  if (data.size() == 0) return;

  if (coords_.record_counts[coords_.n_obs-1]
      + coords_.record_starts[coords_.n_obs-1] > data.size()) {
      std::ostringstream err_stream;
      err_stream << "nemo_feedback::NemoFeedbackWriter::write_variable_profile "
                 << "index range out of bounds '" << variable_name << "' "
                 << coords_.record_counts[coords_.n_obs-1]
                  + coords_.record_starts[coords_.n_obs-1]
                 << " >= " << data.size();
      throw eckit::BadValue(err_stream.str(), Here());
  }
  for (size_t n = 0; n < coords_.n_obs; ++n) {
    var.putVar({n, 0},
               {1, coords_.record_counts[n]},
               data.data()+coords_.record_starts[n]);
  }
}

template <class C>
void NemoFeedbackWriter<C>::write_variable_level_qc(
    const std::string & variable_name,
    const std::vector<int32_t>& data) {
  oops::Log::trace() << "NemoFeedbackWriter::write_variable_level_qc: writing "
                     << variable_name << std::endl;
  auto var = ncFile->getVar(variable_name);
  if (var.isNull()) {
      std::ostringstream err_stream;
      err_stream << "nemo_feedback::NemoFeedbackWriter::write_variable_level_qc"
                 << " ncVar '" << variable_name << "' is not present in "
                 <<"NetCDF file";
      throw eckit::BadValue(err_stream.str(), Here());
    }

  if (data.size() == 0) return;

  if (coords_.record_counts[coords_.n_obs-1]
      + coords_.record_starts[coords_.n_obs-1] > data.size()) {
    std::ostringstream err_stream;
    err_stream << "nemo_feedback::NemoFeedbackWriter::write_variable_level_qc "
               << "index range out of bounds '" << variable_name << "' "
               << coords_.record_counts[coords_.n_obs-1]
                  + coords_.record_starts[coords_.n_obs-1]
               << " >= " << data.size();
    throw eckit::BadValue(err_stream.str(), Here());
  }
  for (size_t n = 0; n < coords_.n_obs; ++n) {
    var.putVar({n, 0},
               {1, coords_.record_counts[n]},
               data.data()+coords_.record_starts[n]);
  }
}

template <class C>
void NemoFeedbackWriter<C>::write_variable_level_qc(
    const std::string & variable_name,
    const std::vector<int32_t>& data,
    const size_t flag_index) {
  oops::Log::trace() << "NemoFeedbackWriter::write_variable_level_qc: writing "
                     << "flag_index: " << flag_index << " of " << variable_name
                     << std::endl;
  auto var = ncFile->getVar(variable_name);
  if (var.isNull()) {
      std::ostringstream err_stream;
      err_stream << "nemo_feedback::NemoFeedbackWriter::write_variable_level_qc"
                 << " ncVar '" << variable_name << "' is not present in "
                 <<"NetCDF file";
      throw eckit::BadValue(err_stream.str(), Here());
    }

  if (data.size() == 0) return;

  if (coords_.record_counts[coords_.n_obs-1]
      + coords_.record_starts[coords_.n_obs-1] > data.size()) {
      std::ostringstream err_stream;
      err_stream << "nemo_feedback::NemoFeedbackWriter::write_variable_level_qc"
                 << " index range out of bounds '" << variable_name << "' "
                 << coords_.record_counts[coords_.n_obs-1]
                  + coords_.record_starts[coords_.n_obs-1]
                 << " >= " << data.size();
      throw eckit::BadValue(err_stream.str(), Here());
  }
  for (size_t n = 0; n < coords_.n_obs; ++n) {
    var.putVar({n, 0, flag_index},
               {1, coords_.record_counts[n], 1},
               data.data()+coords_.record_starts[n]);
  }
}

template class NemoFeedbackWriter<double>;
template class NemoFeedbackWriter<float>;

}  // namespace nemo_feedback
