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

#define QC_CONVENTIONS "U.S. Integrated Ocean Observing System, 2017. Manual "\
  "for the Use of Real-Time Oceanographic Data Quality Control Flags, Version"\
  " 1.1"

namespace nemo_feedback {

const std::map<std::string, size_t> NemoFeedbackWriter::coord_sizes{
  {N_QCF, 2}, {STRINGNAM, STRINGNAM_NUM}, {STRINGGRID, 1},
  {STRINGWMO, STRINGWMO_NUM}, {STRINGTYP, STRINGTYP_NUM}, {STRINGJULD, 14}};

typedef char fixed_length_name_type[STRINGNAM_NUM+1];
typedef char fixed_length_type_type[STRINGTYP_NUM+1];

NemoFeedbackWriter::NemoFeedbackWriter(
    eckit::PathName& filename,
    const size_t n_obs_to_write,
    const std::vector<bool> & to_write,
    const std::vector<double> & lons, 
    const std::vector<double> & lats,
    const std::vector<double> & depths, 
    const std::vector<double> & times,
    const std::vector<std::string> & variable_names,
    const std::vector<std::string> & long_names,
    const std::vector<std::string> & unit_names,
    const std::vector<std::string> & additional_variables,
    const size_t n_levels, 
    const util::DateTime & juld_reference,
    const std::vector<std::string> & station_types )
    : ncFile(nullptr), n_levels_(n_levels) {
  oops::Log::debug() << "nemo_feedback::NemoFieldReader::NemoFieldReader"
                     << " filename: " << filename.fullName().asString()
                     << std::endl;

  ncFile = std::make_unique<netCDF::NcFile>(filename.fullName().asString(),
      netCDF::NcFile::replace);
  if (ncFile->isNull()) {
    std::ostringstream err_stream;
    err_stream << "nemo_feedback::NemoFieldReader::NemoFieldReader cannot open "
               << filename << std::endl;
    eckit::BadValue(err_stream.str(), Here());
  }

  size_t n_obs = lats.size();
  size_t n_obs_vars = variable_names.size();
  size_t max_n_add_entries = additional_variables.size();
  define_coord_variables(
      n_obs_to_write, 
      n_levels, 
      n_obs_vars, 
      max_n_add_entries);
  write_metadata_variables(
      variable_names,  
      additional_variables,
      juld_reference);
  write_coord_variables(
      n_obs,
      n_obs_to_write,
      to_write,
      lons, 
      lats, 
      depths, 
      times);
  define_whole_report_variables();
  write_whole_report_variables(
      n_obs,
      n_obs_to_write,
      to_write, 
      station_types);

  for (int i=0; i < variable_names.size(); ++i) {
    define_variable(
        variable_names[i], 
        long_names[i], 
        unit_names[i],
        additional_variables);
  }
}

void NemoFeedbackWriter::define_coord_variables(
    const size_t n_obs_to_write,
    const size_t n_levels, 
    const size_t n_obs_vars,
    const size_t n_add_entries) {
  netCDF::NcDim tmp_Dim;
  for (const auto& kv : coord_sizes) {
    if (kv.first == N_QCF) {
      nqcf_dim = std::make_unique<netCDF::NcDim>(
          ncFile->addDim(kv.first, kv.second));
    } else {
      tmp_Dim = ncFile->addDim(kv.first, kv.second);
    }
  }

  nobs_dim = std::make_unique<netCDF::NcDim>(ncFile->addDim(N_OBS, n_obs_to_write));
  nlevels_dim = std::make_unique<netCDF::NcDim>(ncFile->addDim(N_LEVELS,
        n_levels));
  tmp_Dim = ncFile->addDim(N_VARS, n_obs_vars);
  tmp_Dim = ncFile->addDim(N_ENTRIES, n_add_entries);
}

void NemoFeedbackWriter::write_metadata_variables(
    const std::vector<std::string>& variable_names,
    const std::vector<std::string>& additional_variables,
    const util::DateTime& juld_reference) {
  {
    netCDF::NcVar nc_juld_var = ncFile->addVar("JULD_REFERENCE", netCDF::ncChar,
        ncFile->getDim(STRINGJULD));
    nc_juld_var.putAtt("long_name", "Date of reference for julian days");
    nc_juld_var.putAtt("Conventions", "YYYYMMDDHHMMSS");
    int year, month, day, hour, minute, second;
    juld_reference.toYYYYMMDDhhmmss(year, month, day, hour, minute, second);
    std::ostringstream ref_stream;
    ref_stream << std::setw(4) << std::setfill('0') << year;
    ref_stream << std::setw(2) << month << day << hour << minute;
    // fudge to handle lack of padding when seconds are 0.
    if (second == 0) {ref_stream << "00";} else {ref_stream << second;};

    nc_juld_var.putVar(ref_stream.str().data());
  }

  {
    size_t n_vars = variable_names.size();
    std::vector<netCDF::NcDim>
        dims{ncFile->getDim(N_VARS), ncFile->getDim(STRINGNAM)};
    netCDF::NcVar nc_var_list_var = ncFile->addVar("VARIABLES",
        netCDF::ncChar, dims);
    nc_var_list_var.putAtt("long_name", "List of variables in feedback files");
    fixed_length_name_type data;
    // trim and pad string to fit in nc char array
    for (size_t i=0; i < n_vars; ++i) {
      for (size_t j=0; j < STRINGNAM_NUM; ++j) {
        if (j < variable_names[i].length()) {
          data[j] = static_cast<char>(variable_names.at(i).at(j));
        } else {data[j] = ' ';}
      }
      nc_var_list_var.putVar({i, 0}, {1, STRINGNAM_NUM}, data);
    }
  }

  {
    size_t max_n_add_entries = additional_variables.size();
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
        if (j < additional_variables[i].length()) {
          data[j] = static_cast<char>(additional_variables.at(i).at(j));
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
    const std::string & units,
    const std::vector<std::string> & additional_names ) {

  const std::vector<netCDF::NcDim> dims{*nobs_dim, *nlevels_dim};

  netCDF::NcVar obs_var = ncFile->addVar(variable_name + "_OBS",
      netCDF::ncDouble, dims);
  obs_var.putAtt("long_name", longName);
  obs_var.putAtt("units", units);

  for (const auto &name : additional_names) {
    netCDF::NcVar add_var = ncFile->addVar(variable_name + "_" + name,
        netCDF::ncDouble, dims);
    add_var.putAtt("long_name", longName + " " + name);
    add_var.putAtt("units", units);
  }

  {
    const std::vector<netCDF::NcDim> qcf_dims{*nobs_dim, *nqcf_dim};
    netCDF::NcVar qc_flags_var = ncFile->addVar(variable_name + "_QC_FLAGS",
        netCDF::ncInt, qcf_dims);
  }

  {
    const std::vector<netCDF::NcDim> lvl_qcf_dims{*nobs_dim, *nlevels_dim,
        *nqcf_dim};
    netCDF::NcVar level_qc_flags_var = ncFile->addVar(variable_name
        + "_LEVEL_QC_FLAGS", netCDF::ncInt, lvl_qcf_dims);
  }

  netCDF::NcVar qc_var = ncFile->addVar(variable_name + "_QC", netCDF::ncInt,
      *nobs_dim);
  qc_var.putAtt("Conventions", QC_CONVENTIONS);

  netCDF::NcVar level_qc_var = ncFile->addVar(variable_name + "_LEVEL_QC",
      netCDF::ncInt, dims);
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

// Routine to subset data before writing.
template <typename T>
std::vector<T> NemoFeedbackWriter::reduce_data(
    const size_t & n_obs,
    const size_t & n_obs_to_write,
    const std::vector<bool> & to_write,
    const std::vector<T> & data_in) {
  std::vector<T> data_out(n_obs_to_write);   
  int j = 0;
  for (int i = 0; i < n_obs; ++i) {
    if (to_write[i]) {
      data_out[j++] = data_in[i];
    }
  }
  return data_out;
}

void NemoFeedbackWriter::write_coord_variables(
    const size_t & n_obs,
    const size_t & n_obs_to_write,
    const std::vector<bool> & to_write,
    const std::vector<double> & lons,
    const std::vector<double> & lats, 
    const std::vector<double> & levels,
    const std::vector<double> & times) {

  netCDF::NcVar lat_var = ncFile->addVar("LATITUDE", netCDF::ncDouble,
      *nobs_dim);
  lat_var.putAtt("units", "degrees_north");
  lat_var.putAtt("long_name", "latitude");
  std::vector<double> reduced_lats = reduce_data(n_obs, 
                                                 n_obs_to_write, 
                                                 to_write, 
                                                 lats);
  lat_var.putVar(reduced_lats.data());

  netCDF::NcVar lon_var = ncFile->addVar("LONGITUDE", netCDF::ncDouble,
      *nobs_dim);
  lon_var.putAtt("units", "degrees_east");
  lon_var.putAtt("long_name", "longitude");
  std::vector<double> reduced_lons = reduce_data(n_obs, 
                                                 n_obs_to_write, 
                                                 to_write, 
                                                 lons);
  lon_var.putVar(reduced_lons.data());

  const std::vector<netCDF::NcDim> dims{*nobs_dim, *nlevels_dim};
  netCDF::NcVar depth_var = ncFile->addVar("DEPTH", netCDF::ncDouble,
      dims);
  depth_var.putAtt("units", "metre");
  depth_var.putAtt("long_name", "Depth");
  std::vector<double> reduced_levels = reduce_data(n_obs, 
                                                   n_obs_to_write, 
                                                   to_write, 
                                                   levels);
  depth_var.putVar(reduced_levels.data());

  netCDF::NcVar juld_var = ncFile->addVar("JULD", netCDF::ncDouble, *nobs_dim);
  juld_var.putAtt("units", "days since JULD_REFERENCE");
  juld_var.putAtt("long_name", "Julian day");
  std::vector<double> reduced_times = reduce_data(n_obs, 
                                                  n_obs_to_write, 
                                                  to_write, 
                                                  times);
  juld_var.putVar(reduced_times.data());
}

void NemoFeedbackWriter::write_whole_report_variables(
    const size_t & n_obs,
    const size_t & n_obs_to_write,
    const std::vector<bool> & to_write,
    const std::vector<std::string> & station_types) {

  // Write station type.    
  auto station_type_var = ncFile->getVar("STATION_TYPE");
  char station_types_char[n_obs_to_write][4];
  size_t nchars = (ncFile->getDim(STRINGTYP)).getSize();
  std::vector<std::string> reduced_station_types = reduce_data(n_obs, 
                                                               n_obs_to_write, 
                                                               to_write, 
                                                               station_types);
  for (int i = 0; i < n_obs_to_write; ++i) {
    strncpy(station_types_char[i], 
            reduced_station_types[i].data(), 
            nchars);
  }
  station_type_var.putVar(station_types_char);

}

void NemoFeedbackWriter::write_variable_surf(
    const size_t & n_obs,
    const size_t & n_obs_to_write,
    const std::vector<bool> & to_write,
    const std::string & variable_name,
    const std::vector<double>& data) {
  oops::Log::trace() << "NemoFeedbackWriter::write_variable_surf: writing "
                     << variable_name << std::endl;
  auto surf_var = ncFile->getVar(variable_name);
  std::vector<double> reduced_data = reduce_data(n_obs, 
                                                 n_obs_to_write, 
                                                 to_write, 
                                                 data);
  surf_var.putVar(reduced_data.data());
}

void NemoFeedbackWriter::write_variable_surf_qc(
    const size_t & n_obs,
    const size_t & n_obs_to_write,
    const std::vector<bool> & to_write,
    const std::string & variable_name, 
    const std::vector<int32_t>& data) {
  oops::Log::trace() << "NemoFeedbackWriter::write_variable_surf_qc: writing "
                     << variable_name << std::endl;
  auto surf_var = ncFile->getVar(variable_name);
  std::vector<int32_t> reduced_data = reduce_data(n_obs, 
                                                  n_obs_to_write, 
                                                  to_write, 
                                                  data);
  surf_var.putVar(reduced_data.data());
}

void NemoFeedbackWriter::write_variable_surf_qc(
    const size_t & n_obs,
    const size_t & n_obs_to_write,
    const std::vector<bool> & to_write,
    const std::string& variable_name, 
    const std::vector<int32_t>& data,
    const size_t flag_index) {
  oops::Log::trace() << "NemoFeedbackWriter::write_variable_surf_qc: writing "
                     << "flag_index: " << flag_index << " of " << variable_name
                     << std::endl;
  auto surf_var = ncFile->getVar(variable_name);
  std::vector<int32_t> reduced_data = reduce_data(n_obs, 
                                                  n_obs_to_write, 
                                                  to_write, 
                                                  data);
  surf_var.putVar({0, flag_index}, {n_obs_to_write, 1}, reduced_data.data());
}

}  // namespace nemo_feedback
