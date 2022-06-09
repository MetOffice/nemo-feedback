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

const std::map<std::string, size_t> NemoFeedbackWriter::coord_sizes{
  {N_QCF, 2}, {STRINGNAM, STRINGNAM_NUM}, {STRINGGRID, 1},
  {STRINGWMO, STRINGWMO_NUM}, {STRINGTYP, STRINGTYP_NUM},
  {STRINGJULD, STRINGJULD_NUM}};

typedef char fixed_length_name_type[STRINGNAM_NUM+1];
typedef char fixed_length_type_type[STRINGTYP_NUM+1];

NemoFeedbackWriter::NemoFeedbackWriter(
    eckit::PathName& filename,
    const size_t & n_obs_to_write,
    const std::vector<bool> & to_write,
    const CoordData & coords,
    const NameData & name_data,
    const std::vector<bool> & extra_vars,
    const std::vector<std::string> & station_types,
    const std::vector<std::string> & station_ids,
    const std::vector<size_t>& record_starts,
    const std::vector<size_t>& record_counts)
    : ncFile(nullptr), coords_(coords), n_obs_(coords.n_obs),
      n_obs_to_write_(n_obs_to_write),
      to_write_(to_write), name_data_(name_data) {
  oops::Log::trace() << "nemo_feedback::NemoFieldReader::NemoFieldReader"
                     << " constructing for: " << filename.fullName().asString()
                     << std::endl;

  ncFile = std::make_unique<netCDF::NcFile>(filename.fullName().asString(),
      netCDF::NcFile::replace);
  if (ncFile->isNull()) {
    std::ostringstream err_stream;
    err_stream << "nemo_feedback::NemoFieldReader::NemoFieldReader cannot open "
               << filename << std::endl;
    eckit::BadValue(err_stream.str(), Here());
  }

  const size_t n_extra = std::count(extra_vars.begin(), extra_vars.end(), true);
  const size_t n_obs_vars = name_data_.variable_names.size() - n_extra;
  const size_t max_n_add_entries = name_data_.additional_names.size();

  if (to_write_.size() != coords.n_locs)
    throw eckit::BadValue(std::string("NemoFeedbackWriter::") +
        " to_write  not defined for some observations to_write_.size() = "
        + std::to_string(to_write_.size()) + " n_locs = "
        + std::to_string(coords.n_locs), Here());

  define_coord_variables(
      n_obs_vars,
      max_n_add_entries,
      n_extra);
  write_metadata_variables(extra_vars);
  write_coord_variables(record_starts, record_counts);
  define_whole_report_variables();
  write_whole_report_variables(
      station_types,
      station_ids);

  for (int i=0; i < name_data_.variable_names.size(); ++i) {
    if (extra_vars[i]) {
      define_extra_variable(name_data_.variable_names[i],
                            name_data_.long_names[i],
                            name_data_.unit_names[i]);
    } else {
      define_variable(
          name_data_.variable_names[i],
          name_data_.long_names[i],
          name_data_.unit_names[i]);
    }
  }
}

void NemoFeedbackWriter::define_coord_variables(
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
      ncFile->addDim(N_OBS, n_obs_to_write_));
  nlevels_dim = std::make_unique<netCDF::NcDim>(ncFile->addDim(N_LEVELS,
        coords_.n_levels));
  tmp_Dim = ncFile->addDim(N_VARS, n_obs_vars);
  tmp_Dim = ncFile->addDim(N_ENTRIES, n_add_entries);
  if (n_extra > 0) {
    tmp_Dim = ncFile->addDim(N_EXTRA, n_extra);
  }
}

void NemoFeedbackWriter::write_metadata_variables(
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
void NemoFeedbackWriter::define_whole_report_variables() {
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
    tmp_var.putAtt("Conventions", "OPS flag conventions");
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
    tmp_var.putAtt("Conventions", "OPS flag conventions");
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
    tmp_var.putAtt("Conventions", "OPS flag conventions");
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
    tmp_var.putAtt("Conventions", "OPS flag conventions");
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
void NemoFeedbackWriter::define_variable(
    const std::string & variable_name,
    const std::string & longName,
    const std::string & units) {

  const std::vector<netCDF::NcDim> dims{*nobs_dim, *nlevels_dim};

  netCDF::NcVar obs_var = ncFile->addVar(variable_name + "_OBS",
      netCDF::ncDouble, dims);
  obs_var.putAtt("long_name", longName);
  obs_var.putAtt("units", units);
  obs_var.setFill(true, double_fillvalue);

  for (const auto &name : name_data_.additional_names) {
    netCDF::NcVar add_var = ncFile->addVar(variable_name + "_" + name,
        netCDF::ncDouble, dims);
    add_var.putAtt("long_name", longName + " " + name);
    add_var.putAtt("units", units);
    add_var.setFill(true, double_fillvalue);
  }

  {
    const std::vector<netCDF::NcDim> qcf_dims{*nobs_dim, *nqcf_dim};
    netCDF::NcVar qc_flags_var = ncFile->addVar(variable_name + "_QC_FLAGS",
        netCDF::ncInt, qcf_dims);
    qc_flags_var.setFill(true, int32_fillvalue);
    qc_flags_var.putAtt("long_name", std::string("quality flags on ")
                        + longName);
    qc_flags_var.putAtt("Conventions", "OPS flag conventions");
  }

  {
    const std::vector<netCDF::NcDim> lvl_qcf_dims{*nobs_dim, *nlevels_dim,
        *nqcf_dim};
    netCDF::NcVar level_qc_flags_var = ncFile->addVar(variable_name
        + "_LEVEL_QC_FLAGS", netCDF::ncInt, lvl_qcf_dims);
    level_qc_flags_var.setFill(true, int32_fillvalue);
    level_qc_flags_var.putAtt("long_name",
        std::string("quality flags for each level on ") + longName);
    level_qc_flags_var.putAtt("Conventions", "OPS flag conventions");
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
void NemoFeedbackWriter::define_extra_variable(
    const std::string & variable_name,
    const std::string & longName,
    const std::string & units ) {

  const std::vector<netCDF::NcDim> dims{*nobs_dim, *nlevels_dim};

  netCDF::NcVar var = ncFile->addVar(variable_name,
                                     netCDF::ncDouble, dims);
  var.putAtt("long_name", longName);
  var.putAtt("units", units);
}

template <typename T>
std::vector<T> NemoFeedbackWriter::reduce_data(
    const std::vector<T> & data_in) {
  // profile data, n_obs dimension var - using the n_obs should be fine
  // surface data, n_locs/n_obs dimension var, n_locs = n_obs
  // so using n_obs should be fine
  std::vector<T> data_out(n_obs_to_write_);
  int j = 0;
  for (int i = 0; i < n_obs_; ++i) {
    if (to_write_[i]) {
      data_out[j++] = data_in[i];
    }
  }
  return data_out;
}

template <typename T>
void NemoFeedbackWriter::reduce_profile_data(
    const std::vector<size_t> & record_starts,
    const std::vector<size_t> & record_counts,
    const std::vector<T> & data_in,
    std::vector<size_t> & record_starts_out,
    std::vector<size_t> & record_counts_out,
    std::vector<T> & data_out
    ) {
  // with profile data n_obs != n_locs, and so we setup new record_starts and
  // counts based on the new data vector.
  data_out.reserve((std::count(to_write_.begin(), to_write_.end(), true)));
  record_starts_out.reserve(n_obs_);
  record_counts_out.reserve(n_obs_);
  for (int i = 0; i < n_obs_; ++i) {
    size_t reclen = 0;
    for (int l = record_starts[i]; l < record_starts[i]+record_counts[i]; ++l) {
      if (to_write_[l]) {
        if (i >= record_starts_out.size()) {
          if (i == 0) {
            record_starts_out.push_back(l);
          } else {
            record_starts_out.push_back(
                record_starts_out[i-1]+record_counts_out[i-1]);
          }
        }
        data_out.push_back(data_in[l]);
        reclen++;
      }
    }
    record_counts_out.push_back(reclen);
  }
}

void NemoFeedbackWriter::write_coord_variables(
    const std::vector<size_t>& record_starts,
    const std::vector<size_t>& record_counts) {

  netCDF::NcVar lat_var = ncFile->addVar("LATITUDE", netCDF::ncDouble,
      *nobs_dim);
  lat_var.putAtt("units", "degrees_north");
  lat_var.putAtt("long_name", "latitude");
  if (n_obs_ == coords_.n_locs) {
    std::vector<double> reduced_lats = reduce_data(coords_.lats);
    lat_var.putVar(reduced_lats.data());
  } else {
    lat_var.putVar(coords_.lats.data());
  }

  netCDF::NcVar lon_var = ncFile->addVar("LONGITUDE", netCDF::ncDouble,
      *nobs_dim);
  lon_var.putAtt("units", "degrees_east");
  lon_var.putAtt("long_name", "longitude");
  if (n_obs_ == coords_.n_locs) {
      std::vector<double> reduced_lons = reduce_data(coords_.lons);
      lon_var.putVar(reduced_lons.data());
  } else {
      lon_var.putVar(coords_.lons.data());
  }

  const std::vector<netCDF::NcDim> dims{*nobs_dim, *nlevels_dim};
  netCDF::NcVar depth_var = ncFile->addVar("DEPTH", netCDF::ncDouble,
      dims);
  depth_var.putAtt("units", "metre");
  depth_var.putAtt("long_name", "Depth");
  if (n_obs_ == coords_.n_locs) {
    std::vector<double> reduced_depths = reduce_data(coords_.depths);
    depth_var.putVar(reduced_depths.data());
  } else {
    std::vector<size_t> reduced_record_starts, reduced_record_counts;
    std::vector<double> reduced_depths;
    reduce_profile_data(record_starts, record_counts, coords_.depths,
        reduced_record_starts, reduced_record_counts, reduced_depths);
    for (size_t n = 0; n < n_obs_; ++n) {
      depth_var.putVar({n, 0},
                       {1, reduced_record_counts[n]},
                       reduced_depths.data()+reduced_record_starts[n]);
    }
  }

  netCDF::NcVar juld_var = ncFile->addVar("JULD", netCDF::ncDouble, *nobs_dim);
  juld_var.putAtt("units", "days since JULD_REFERENCE");
  juld_var.putAtt("long_name", "Julian day");
  if (n_obs_ == coords_.n_locs) {
    std::vector<double> reduced_times = reduce_data(coords_.julian_days);
    juld_var.putVar(reduced_times.data());
  } else {
    juld_var.putVar(coords_.julian_days.data());
  }
}

void NemoFeedbackWriter::write_whole_report_variables(
    const std::vector<std::string> & station_types,
    const std::vector<std::string> & station_ids) {

  if (station_ids.size() != n_obs_ ||
      station_types.size() != n_obs_ )
    throw eckit::BadValue(std::string("NemoFeedbackWriter::") +
        "write_whole_report_variables: station_ids or station_types " +
        "not defined for some observations", Here());

  // Write station type.
  {
    size_t nchars = (ncFile->getDim(STRINGTYP)).getSize();
    auto station_type_var = ncFile->getVar("STATION_TYPE");
    int j = 0;
    // +1 is for the null-terminator of a cstring
    char* buffer = new char[n_obs_to_write_*nchars+1];
    for (int i = 0; i < n_obs_; ++i) {
      if (n_obs_ == coords_.n_locs) {
        if (to_write_[i]) strcpy(buffer + nchars*j++, station_types[i].c_str());
      } else {
        strcpy(buffer + nchars*j++, station_types[i].c_str());
      }
    }
    station_type_var.putVar({0, 0},
                            {n_obs_to_write_, nchars},
                            buffer);
    delete[] buffer;
  }

  // Write station IDs.
  {
    size_t nchars = (ncFile->getDim(STRINGWMO)).getSize();
    auto station_id_var = ncFile->getVar("STATION_IDENTIFIER");
    int j = 0;
    // +1 is for the null-terminator of a cstring
    char* buffer = new char[n_obs_to_write_*nchars+1];
    for (int i = 0; i < n_obs_; ++i) {
      if (n_obs_ == coords_.n_locs) {
        if (to_write_[i]) strcpy(buffer + nchars*j++, station_ids[i].c_str());
      } else {
        strcpy(buffer + nchars*j++, station_ids[i].c_str());
      }
    }
    station_id_var.putVar({0, 0},
                          {n_obs_to_write_, nchars},
                          buffer);
    delete[] buffer;
  }
}

void NemoFeedbackWriter::write_variable_surf(
    const std::string & variable_name,
    const std::vector<double>& data) {
  oops::Log::trace() << "NemoFeedbackWriter::write_variable_surf: writing "
                     << variable_name << std::endl;
  auto surf_var = ncFile->getVar(variable_name);
  std::vector<double> reduced_data = reduce_data(data);
  surf_var.putVar(reduced_data.data());
}

void NemoFeedbackWriter::write_variable_surf_qc(
    const std::string & variable_name,
    const std::vector<int32_t>& data) {
  oops::Log::trace() << "NemoFeedbackWriter::write_variable_surf_qc: writing "
                     << variable_name << std::endl;
  auto surf_var = ncFile->getVar(variable_name);
  std::vector<int32_t> reduced_data = reduce_data(data);
  surf_var.putVar(reduced_data.data());
}

void NemoFeedbackWriter::write_variable_surf_qc(
    const std::string& variable_name,
    const std::vector<int32_t>& data,
    const size_t flag_index) {
  oops::Log::trace() << "NemoFeedbackWriter::write_variable_surf_qc: writing "
                     << "flag_index: " << flag_index << " of " << variable_name
                     << std::endl;
  auto surf_var = ncFile->getVar(variable_name);
  std::vector<int32_t> reduced_data = reduce_data(data);
  surf_var.putVar({0, flag_index}, {n_obs_to_write_, 1}, reduced_data.data());
}

void NemoFeedbackWriter::write_variable_profile(
    const std::string & variable_name,
    const std::vector<double>& data,
    const std::vector<size_t>& record_starts,
    const std::vector<size_t>& record_counts) {
  oops::Log::trace() << "NemoFeedbackWriter::write_variable_profile: writing "
                     << variable_name << std::endl;
  auto var = ncFile->getVar(variable_name);
  if (var.isNull()) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFeedbackWriter::write_variable_profile "
                 << "ncVar '" << variable_name << "' is not present in "
                 <<"NetCDF file";
      throw eckit::BadValue(err_stream.str(), Here());
    }

  std::vector<size_t> reduced_record_starts, reduced_record_counts;
  std::vector<double> reduced_data;
  reduce_profile_data(record_starts, record_counts, data, reduced_record_starts,
      reduced_record_counts, reduced_data);

  for (size_t n = 0; n < n_obs_; ++n) {
    var.putVar({n, 0},
               {1, reduced_record_counts[n]},
               reduced_data.data()+reduced_record_starts[n]);
  }
}

void NemoFeedbackWriter::write_variable_level_qc(
    const std::string & variable_name,
    const std::vector<int32_t>& data,
    const std::vector<size_t>& record_starts,
    const std::vector<size_t>& record_counts) {
  oops::Log::trace() << "NemoFeedbackWriter::write_variable_level_qc: writing "
                     << variable_name << std::endl;
  auto var = ncFile->getVar(variable_name);
  if (var.isNull()) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFeedbackWriter::write_variable_level_qc "
                 << "ncVar '" << variable_name << "' is not present in "
                 <<"NetCDF file";
      throw eckit::BadValue(err_stream.str(), Here());
    }

  std::vector<size_t> reduced_record_starts, reduced_record_counts;
  std::vector<int32_t> reduced_data;
  reduce_profile_data(record_starts, record_counts, data, reduced_record_starts,
      reduced_record_counts, reduced_data);

  for (size_t n = 0; n < n_obs_; ++n) {
    var.putVar({n, 0},
               {1, reduced_record_counts[n]},
               reduced_data.data()+reduced_record_starts[n]);
  }
}

void NemoFeedbackWriter::write_variable_level_qc(
    const std::string & variable_name,
    const std::vector<int32_t>& data,
    const size_t flag_index,
    const std::vector<size_t>& record_starts,
    const std::vector<size_t>& record_counts) {
  oops::Log::trace() << "NemoFeedbackWriter::write_variable_level_qc: writing "
                     << "flag_index: " << flag_index << " of " << variable_name
                     << std::endl;
  auto var = ncFile->getVar(variable_name);
  if (var.isNull()) {
      std::ostringstream err_stream;
      err_stream << "orcamodel::NemoFeedbackWriter::write_variable_level_qc "
                 << "ncVar '" << variable_name << "' is not present in "
                 <<"NetCDF file";
      throw eckit::BadValue(err_stream.str(), Here());
    }

  std::vector<size_t> reduced_record_starts, reduced_record_counts;
  std::vector<int32_t> reduced_data;
  reduce_profile_data(record_starts, record_counts, data, reduced_record_starts,
      reduced_record_counts, reduced_data);

  for (size_t n = 0; n < n_obs_; ++n) {
    var.putVar({n, 0, flag_index},
               {1, reduced_record_counts[n], 1},
               reduced_data.data()+reduced_record_starts[n]);
  }
}

}  // namespace nemo_feedback
