/*
 * (C) Crown Copyright 2020, the Met Office. All rights reserved.
 *
 * Refer to COPYRIGHT.txt of this distribution for details.
 */

#include "nemo-feedback/NemoFeedback.h"

#include <string.h>

#include <utility>
#include <vector>
#include <bitset>
#include <set>

#include "eckit/mpi/Comm.h"
#include "eckit/mpi/Parallel.h"
#include "eckit/exception/Exceptions.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/util/Logger.h"
#include "oops/util/DateTime.h"
#include "oops/util/missingValues.h"
#include "nemo-feedback/NemoFeedbackParameters.h"
#include "nemo-feedback/NemoFeedbackWriter.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/filters/DiagnosticFlag.h"
#include "ufo/filters/processWhere.h"
#include "ufo/utils/metoffice/MetOfficeQCFlags.h"


namespace nemo_feedback {

NemoFeedback::NemoFeedback(
    ioda::ObsSpace & obsdb,
    const Parameters_ & params,
    std::shared_ptr< ioda::ObsDataVector<int> > flags,
    std::shared_ptr< ioda::ObsDataVector<float> > obsErrors)
      :
    obsdb_(obsdb),
    data_(obsdb_),
    geovars_(),
    flags_(std::move(flags)),
    obsErrors_(std::move(obsErrors)),
    parameters_(params),
    validityTime_(obsdb.windowStart() +
        (obsdb.windowEnd() - obsdb.windowStart()) / 2)
{
  oops::Log::trace() << "NemoFeedback constructor starting" << std::endl;

  const std::vector<int> channels{};
  std::vector<std::string> varnames;
  for (const NemoFeedbackVariableParameters& nemoVariableParams :
        parameters_.variables.value()) {
    const std::string varname = nemoVariableParams.name.value();
    if (std::find(varnames.begin(), varnames.end(), varname) == varnames.end())
      varnames.push_back(varname);
  }
  geovars_ = oops::Variables(varnames, channels);

  if (obsdb.comm().size() > 1)
    throw eckit::BadValue("nemo-feedback filter not MPI capable", Here());
  MPI_Comm mpiComm = MPI_COMM_WORLD;
  if (auto parallelComm =
        dynamic_cast<const eckit::mpi::Parallel*>(&obsdb.comm())) {
    mpiComm = parallelComm->MPIComm();
  }
}

NemoFeedback::~NemoFeedback() {
  oops::Log::trace() << "NemoFeedback destructor " << std::endl;
}

void NemoFeedback::priorFilter(const ufo::GeoVaLs & gv) {
  oops::Log::trace() << "NemoFeedback priorFilter" << std::endl;
}

void NemoFeedback::postFilter(const ufo::GeoVaLs & gv,
                  const ioda::ObsVector &ov,
                  const ioda::ObsVector &bv,
                  const ufo::ObsDiagnostics &dv) {
  oops::Log::trace() << "NemoFeedback postFilter" << std::endl;

  eckit::PathName test_data_path(parameters_.Filename);

  // Generate lists of the variable names to setup the file
  std::vector<std::string> additional_names;
  std::vector<std::string> variable_names;
  std::vector<std::string> long_names;
  std::vector<std::string> unit_names;
  std::vector<bool> extra_vars;
  bool is_altimeter = false;
  bool is_profile = false;
  for (const NemoFeedbackVariableParameters& nemoVariableParams :
        parameters_.variables.value()) {
    if (nemoVariableParams.nemoName.value() == "SLA") {
      is_altimeter = true;
    }
    if (nemoVariableParams.nemoName.value() == "POTM" ||
        nemoVariableParams.nemoName.value() == "PSAL") {
      is_profile = true;
    }
    variable_names.push_back(nemoVariableParams.nemoName.value());
    long_names.push_back(nemoVariableParams.longName.value());
    unit_names.push_back(nemoVariableParams.units.value());
    extra_vars.push_back(nemoVariableParams.extravar.value().value_or(false));
    auto additionalVariablesParams = nemoVariableParams.variables.value();
    for (const NemoFeedbackAddVariableParameters& addVariableParams :
        additionalVariablesParams) {
      auto add_suffix = addVariableParams.feedbackSuffix.value();
      if (std::find(additional_names.begin(), additional_names.end(),
            add_suffix) == additional_names.end()) {
        additional_names.push_back(add_suffix);
      }
    }
  }

  if (is_profile && is_altimeter)
    throw eckit::BadValue(std::string("NemoFeedback::postFilter cannot write")
        + " profile and altimeter data to the same file", Here());

  //  Calculate n_obs, n_levels, starts, counts
  //  Fill out lats, lons, depths, julian_days
  // ov.n_obs is the number of non-missing non-qced obs?
  // obsdb_.nlocs is the number of original point-observation locations.
  // TODO: We need lats/lons/depths to be of size N_OBS for profiles

  NemoFeedbackWriter::CoordData coords;
  std::vector<size_t> record_starts;
  std::vector<size_t> record_counts;
  coords.n_levels = 1;

  groupCoordsByRecord(obsdb_.nlocs(), coords, record_starts,
      record_counts, is_profile);

  size_t n_locs = obsdb_.nlocs();
  size_t n_obs = obsdb_.nlocs();
  if (is_profile) {
    n_obs = record_starts.size();
  }

  // Handle the where option.
  std::vector<bool> to_write = ufo::processWhere(parameters_.where, data_,
                                   ufo::WhereOperator::AND);
  if (to_write.size() != n_obs) {
    to_write.assign(n_obs, true);
  }
  // Calculate total number of obs to actually write.
  const size_t n_obs_to_write = std::count(to_write.begin(), to_write.end(),
      true);


  // Set up the station type and identifier variables. The way to do this
  // depends on the data type. Also modify which obs to write if this is
  // altimeter data.
  std::vector<std::string> station_types(n_obs, "    ");
  std::vector<std::string> station_ids(n_obs, "        ");

  if (is_altimeter) {
    setupAltimeterIds(n_obs, station_ids, station_types, to_write);
  } else {
    setupIds(n_obs, record_starts, record_counts, station_ids, station_types);
  }

  NemoFeedbackWriter fdbk_writer(
      test_data_path,
      n_obs_to_write,
      to_write,
      coords.lons,
      coords.lats,
      coords.depths,
      coords.julian_days,
      variable_names,
      long_names,
      unit_names,
      additional_names,
      extra_vars,
      coords.n_levels,
      coords.juld_reference,
      station_types,
      station_ids);

  // Write the data
  std::vector<double> variable_data;
  std::vector<int> variable_qcFlags(n_obs, 0);
  std::vector<int> variable_qc(n_obs);
  std::vector<int> variable_level_qc(n_locs);
  std::vector<ufo::DiagnosticFlag> final_qc;
  std::vector<ufo::DiagnosticFlag> do_not_assimilate;
  for (const NemoFeedbackVariableParameters& nemoVariableParams :
        parameters_.variables.value()) {
    auto nemo_name = nemoVariableParams.nemoName.value();
    auto ufo_name = nemoVariableParams.name.value();
    obsdb_.get_db("ObsValue", ufo_name, variable_data);
    auto missing_value = util::missingValue(variable_data[0]);
    oops::Log::trace() << "Missing value OBS: " << missing_value << std::endl;
    for (int i=0; i < n_locs; ++i) {
      if (variable_data[i] == missing_value)
        variable_data[i] = NemoFeedbackWriter::double_fillvalue;
    }
    auto extra_var = nemoVariableParams.extravar.value().value_or(false);
    if (extra_var) {
      fdbk_writer.write_variable_surf(
          n_obs,
          n_obs_to_write,
          to_write,
          nemo_name,
          variable_data);
      // If this is an extra variable we do not want to write any of the
      // other variables with _OBS, _QC etc. added on to the name.
      continue;
    }

    if (is_profile) {
      fdbk_writer.write_variable_profile(
          n_obs,
          nemo_name + "_OBS",
          variable_data,
          record_starts,
          record_counts);
    } else {
      fdbk_writer.write_variable_surf(
          n_obs,
          n_obs_to_write,
          to_write,
          nemo_name + "_OBS",
          variable_data);
    }

    // Write Met Office QC flag data for this variable if they exist.
    if (obsdb_.has("QCFlags", ufo_name)) {
      obsdb_.get_db("QCFlags", ufo_name, variable_qcFlags);

      // To the first qc flag index location
      if (is_profile) {
        fdbk_writer.write_variable_level_qc(
            n_obs,
            nemo_name + "_LEVEL_QC_FLAGS",
            variable_qcFlags,
            0,
            record_starts,
            record_counts);
      } else {
        fdbk_writer.write_variable_surf_qc(
            n_obs,
            n_obs_to_write,
            to_write,
            nemo_name + "_QC_FLAGS",
            variable_qcFlags, 0);
      }
    }

    // Whole Observation report QC flags
    for (int i=0; i < variable_qc.size(); ++i) {
      if (variable_qcFlags[i]
          & ufo::MetOfficeQCFlags::WholeObReport::FinalRejectReport) {
        variable_qc[i] = 4;
      } else {variable_qc[i] = 0;}
    }
    fdbk_writer.write_variable_surf_qc(
        n_obs,
        n_obs_to_write,
        to_write,
        "OBSERVATION_QC",
        variable_qc);

    // Overall quality control flags
    obsdb_.get_db("DiagnosticFlags/FinalReject", ufo_name, final_qc);
    for (int i=0; i < n_locs; ++i) {
      if (final_qc[i]) {
        variable_level_qc[i] = 4;
      } else {variable_level_qc[i] = 1;}
    }
    for (int i=0; i < n_obs; ++i) {
      int j = record_starts[i];
      if (final_qc[j]) {
        variable_qc[j] = 4;
      } else {variable_qc[j] = 1;}
    }
    // Add do not assimilate flag if required.
    if (obsdb_.has("DiagnosticFlags/DoNotAssimilate", ufo_name)) {
      obsdb_.get_db("DiagnosticFlags/DoNotAssimilate", ufo_name,
          do_not_assimilate);
      for (int i=0; i < n_locs; ++i) {
        if (do_not_assimilate[i]) {
          variable_level_qc[i] += 128;
        }
      }
      for (int i=0; i < n_obs; ++i) {
        int j = record_starts[i];
        if (do_not_assimilate[j]) {
          variable_qc[j] += 128;
        }
      }
    }
    fdbk_writer.write_variable_surf_qc(
        n_obs,
        n_obs_to_write,
        to_write,
        nemo_name + "_QC",
        variable_qc);
    if (is_profile) {
      fdbk_writer.write_variable_level_qc(
          n_obs,
          nemo_name + "_LEVEL_QC",
          variable_level_qc,
          record_starts,
          record_counts);
    } else {
      fdbk_writer.write_variable_surf_qc(
          n_obs,
          n_obs_to_write,
          to_write,
          nemo_name + "_LEVEL_QC",
          variable_level_qc);
    }

    // Write additional variables for this variable
    auto additionalVariablesParams = nemoVariableParams.variables.value();
    for (const NemoFeedbackAddVariableParameters& addParams :
        additionalVariablesParams) {
      auto add_name = nemo_name + "_" + addParams.feedbackSuffix.value();
      std::string ioda_group = addParams.iodaGroup.value();
      if (ioda_group == "HofX") {
        oops::Log::trace() << ov << std::endl;
        if (ov.has(ufo_name)) {
          std::vector<std::string> ov_varnames = ov.varnames().variables();
          auto var_it = std::find(ov_varnames.begin(), ov_varnames.end(),
              ufo_name);
          oops::Log::trace() << "NemoFeedback::ov_varnames = "
                             << ov_varnames << std::endl;
          std::size_t var_it_dist = static_cast<std::size_t>(
                                    std::distance(ov_varnames.begin(), var_it));
          oops::Log::trace() << "NemoFeedback::iterator distance is "
                             << var_it_dist
                             << " ov.nvars: " << ov.nvars()
                             << std::endl;
          const auto missing_value_add = util::missingValue(ov[var_it_dist]);
          oops::Log::trace() << "NemoFeedback::Missing value: "
                             << missing_value_add << std::endl;
          for (int i=0; i < n_obs; ++i) {
            const size_t indx = i * ov.nvars() + var_it_dist;
            if (ov[indx] == missing_value_add) {
              variable_data[i] = NemoFeedbackWriter::double_fillvalue;
            } else {
              variable_data[i] = ov[indx];
            }
          }
          if (is_profile) {
            fdbk_writer.write_variable_profile(
                n_obs,
                add_name,
                variable_data,
                record_starts,
                record_counts);
          } else {
            fdbk_writer.write_variable_surf(
                n_obs,
                n_obs_to_write,
                to_write,
                add_name,
                variable_data);
          }
        }
      } else if (!obsdb_.has(ioda_group, ufo_name)) {
        oops::Log::trace() << obsdb_ << std::endl;
        throw eckit::BadValue("NemoFeedback::postFilter missing variable: "
            + ioda_group + "/" + ufo_name, Here());
      } else {
        obsdb_.get_db(ioda_group, ufo_name, variable_data);
        auto missing_value = util::missingValue(variable_data[0]);
        for (int i=0; i < n_obs; ++i) {
          if (variable_data[i] == missing_value)
            variable_data[i] = NemoFeedbackWriter::double_fillvalue;
        }
        if (is_profile) {
          fdbk_writer.write_variable_profile(
                n_obs,
                add_name,
                variable_data,
                record_starts,
                record_counts);
        } else {
          fdbk_writer.write_variable_surf(
              n_obs,
              n_obs_to_write,
              to_write,
              add_name,
              variable_data);
        }
      }
    }
  }
}

void NemoFeedback::groupCoordsByRecord(size_t n_obs,
                                       NemoFeedbackWriter::CoordData& coords,
                                       std::vector<size_t>& record_starts,
                                       std::vector<size_t>& record_counts,
                                       bool is_profile) const {
    size_t n_locs = obsdb_.nlocs();
    if (n_obs != n_locs)
      throw eckit::BadValue(std::string("NemoFeedback::setupJulianDay: ")
          + "n_obs != nlocs (doesn't support observations at depth)",
          Here());

    std::vector<util::DateTime> datetimes(n_locs);
    obsdb_.get_db("MetaData", "dateTime", datetimes);
    coords.lats.resize(n_locs);
    obsdb_.get_db("MetaData", "latitude", coords.lats);
    coords.lons.resize(n_locs);
    obsdb_.get_db("MetaData", "longitude", coords.lons);

    // Handle the reference date option.
    coords.juld_reference = parameters_.refDate.value().value_or(datetimes[0]);

    if (is_profile) {
      std::vector<double> record_lats(1, coords.lats[0]);
      std::vector<double> record_lons(1, coords.lons[0]);
      std::vector<util::DateTime> record_dts(1, datetimes[0]);
      record_counts.resize(1);
      record_starts.resize(1);
      int r = 0;
      n_obs = 1;
      for (int l=0; l<n_locs; ++l) {
        if ((std::abs(coords.lats[l] - record_lats[r]) > 1e-6) ||
            (std::abs(coords.lons[l] - record_lons[r]) > 1e-6) ||
            (std::abs((datetimes[l] - record_dts[r]).toSeconds()) > 0)) {
          record_lats.push_back(coords.lats[l]);
          record_lons.push_back(coords.lons[l]);
          record_dts.push_back(datetimes[l]);
          record_counts.push_back(1);
          record_starts.push_back(l);
          r++;
          n_obs++;
        } else {
          record_counts[r]++;
        }
      }
      // n_levels is equal to largest number of levels across an entire
      // profile data.
      coords.n_levels = *std::max_element(record_counts.begin(),
          record_counts.end());
      coords.lats.resize(n_obs);
      coords.lats.assign(record_lats.begin(), record_lats.end());
      coords.lons.resize(n_obs);
      coords.lons.assign(record_lons.begin(), record_lons.end());
      datetimes.resize(n_obs);
      datetimes.assign(record_dts.begin(), record_dts.end());
      coords.depths.resize(n_locs);
      obsdb_.get_db(parameters_.depthGroup.value().value_or("MetaData"),
                    parameters_.depthVariable.value().value_or("depth_m"),
                    coords.depths);
    } else {
      coords.depths.resize(n_obs);
      std::fill(coords.depths.begin(), coords.depths.end(), 0);
      record_counts.assign(n_obs, 1);
      record_starts.resize(n_obs);
      for (int l=0; l<n_locs; ++l)
        record_starts[l] = l;
    }

    coords.julian_days.resize(n_locs);
    for (int i=0; i < n_obs; ++i) {
      util::Duration duration = static_cast<util::DateTime>(datetimes[i])
          - coords.juld_reference;
      coords.julian_days[i] = duration.toSeconds() / 86400.0;
    }
}

void NemoFeedback::setupIds(const size_t n_obs,
    const std::vector<size_t>& record_starts,
    const std::vector<size_t>& record_counts,
    std::vector<std::string>& station_ids,
    std::vector<std::string>& station_types) const {
    //if (n_obs != obsdb_.nlocs())
    //  throw eckit::BadValue(std::string("NemoFeedback::setupIds ")
    //      + "n_obs != nlocs. Surface observations don't have depth",
    //      Here());
    // Define the station type variable. This may not always be defined in
    // obsdb_.
    if (obsdb_.has("MetaData", "fdbk_station_type")) {
      std::vector<int> station_types_int(obsdb_.nlocs());
      obsdb_.get_db("MetaData", "fdbk_station_type", station_types_int);
      char buffer[5];
      for (int i=0; i < n_obs; ++i) {
        int j = record_starts[i];
        snprintf(buffer, sizeof(buffer), "%-4d", station_types_int[j]);
        station_types[i] = buffer;
      }
    }

    // Define the station identifier variable. These may not always be defined
    // in odbsb_ and may be contained in either the station_id or
    // buoy_identifier variables, or both.
    if (obsdb_.has("MetaData", "station_id")) {
      std::vector<std::string> station_ids_tmp(obsdb_.nlocs());
      obsdb_.get_db("MetaData", "station_id", station_ids_tmp);
      for (int i=0; i< n_obs; ++i) {
        if (station_ids_tmp[i] != "") {
               int j = record_starts[i];
               station_ids_tmp[j].resize(8, ' ');
               station_ids[i] = station_ids_tmp[j].substr(0, 8);
        }
      }
    }
    if (obsdb_.has("MetaData", "buoy_identifier")) {
      std::vector<int> buoy_ids(obsdb_.nlocs());
      obsdb_.get_db("MetaData", "buoy_identifier", buoy_ids);
      char buffer[9];
      auto buoy_id_missing_value = util::missingValue(buoy_ids[0]);
      for (int i=0; i < n_obs; ++i) {
        if (buoy_ids[i] != buoy_id_missing_value) {
          int j = record_starts[i];
          snprintf(buffer, sizeof(buffer), "%-8d", buoy_ids[j]);
          station_ids[i] = buffer;
        }
      }
    }
}

void NemoFeedback::setupAltimeterIds(const size_t n_obs,
    std::vector<std::string>& station_ids,
    std::vector<std::string>& station_types,
    std::vector<bool>& to_write) const {
    if (n_obs != obsdb_.nlocs())
      throw eckit::BadValue(std::string("NemoFeedback::setupAltimeterIds ")
          + "n_obs != nlocs. Altimeter observations don't have depth",
          Here());

    std::vector<util::DateTime> datetimes(n_obs);
    obsdb_.get_db("MetaData", "dateTime", datetimes);

    // Station type and station identifier are both defined from the
    // satellite_identifier metadata.
    std::vector<int> satellite_ids(n_obs);
    obsdb_.get_db("MetaData", "satellite_identifier", satellite_ids);

    // Get the most recent versions of the data.
    // Read version information.
    std::vector<int> version(n_obs);
    obsdb_.get_db("MetaData", "instrument_type", version);
    char buffer9[9];
    char buffer5[5];
    for (int i=0; i < n_obs; ++i) {
      snprintf(buffer9, sizeof(buffer9), "%04d    ", satellite_ids[i]);
      station_ids[i] = buffer9;
      snprintf(buffer5, sizeof(buffer5), "%4d", satellite_ids[i]);
      station_types[i] = buffer5;
    }

    // Get integer values for the day of each data point.
    std::vector<int> ymd(n_obs);
    int ymd_tmp;
    int hms_tmp;
    for (int i=0; i < n_obs; ++i) {
      datetimes[i].toYYYYMMDDhhmmss(ymd_tmp, hms_tmp);
      ymd[i] = ymd_tmp;
    }

    // Get unique values.
    std::set<int> ymd_set(ymd.begin(), ymd.end());
    std::set<int> sid_set(satellite_ids.begin(), satellite_ids.end());

    // Loop through the unique values of ymd and satellite_ids.
    int latest_version;
    auto version_missing_value = util::missingValue(version[0]);
    for (int ymd_to_find : ymd_set) {
      for (int sid_to_find : sid_set) {
        latest_version = 0;
        for (int i=0; i < n_obs; ++i) {
          if (to_write[i] &&
              ymd[i] == ymd_to_find &&
              satellite_ids[i] == sid_to_find &&
              version[i] != version_missing_value &&
              version[i] > latest_version) {
            latest_version = version[i];
          }
        }
        // If latest_version is still 0, this means that either there are
        // no data with a version number defined or all the data are marked
        // as not to write already. Therefore, we only need to do anything
        // further if latest_version > 0.
        if (latest_version > 0) {
          for (int i=0; i < n_obs; ++i) {
            if (to_write[i] &&
                ymd[i] == ymd_to_find &&
                satellite_ids[i] == sid_to_find &&
                version[i] < latest_version) {
              to_write[i] = false;
            }
          }
        }
      }
    }
}

void NemoFeedback::print(std::ostream & os) const {
  os << "NemoFeedback: config = " << parameters_ << std::endl;
}

}  // namespace nemo_feedback
