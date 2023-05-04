/*
 * (C) Crown Copyright 2020, the Met Office. All rights reserved.
 *
 * Refer to COPYRIGHT.txt of this distribution for details.
 */

#include "nemo-feedback/NemoFeedback.h"

#include <string.h>

#include <algorithm>
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
#include "nemo-feedback/NemoFeedbackReduce.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/filters/DiagnosticFlag.h"
#include "ufo/filters/processWhere.h"
#include "ufo/filters/ObsAccessor.h"
#include "ufo/utils/metoffice/MetOfficeQCFlags.h"


namespace nemo_feedback {

/// \brief interface between JEDI and NemoFeedbackWriter to manage the reading
///        and writing of observation and HofX data.
template<typename T>
struct VariableData {
  std::vector<T> data;
  std::vector<T> reduced;

  void read_obs(const ioda::ObsSpace & obsdb, const std::string& obs_group,
      const std::string& ufo_name) {
      obsdb.get_db(obs_group, ufo_name, data);
      T missing_value = util::missingValue(T(0));
      oops::Log::trace() << "Missing value OBS: " << missing_value
                         << std::endl;
      // Convert oops missing values to NEMO missing value
      for_each(data.begin(), data.end(),
          [missing_value](T& d){
            if (d == missing_value) {d = typeToFill::value<T>();} });
  }
  void read_hofx(const ioda::ObsVector &ov, const CoordData &coords,
      const std::string& ufo_name) {
      std::vector<std::string> ov_varnames = ov.varnames().variables();
      auto var_it = std::find(ov_varnames.begin(), ov_varnames.end(),
          ufo_name);
      oops::Log::trace() << "NemoFeedback::variableData begin ov_varnames = "
                         << ov_varnames << std::endl;
      std::size_t var_it_dist = static_cast<std::size_t>(
                                std::distance(ov_varnames.begin(), var_it));
      oops::Log::debug() << "NemoFeedback::variableData::iterator distance is "
                         << var_it_dist << " ov.nvars: " << ov.nvars()
                         << std::endl;
      // Convert oops missing values to NEMO missing value
      double missing_value = util::missingValue(static_cast<double>(0));
      if (coords.n_locs > 0)
          missing_value = util::missingValue(ov[var_it_dist]);
      oops::Log::trace() << "Missing value HofX: " << missing_value
                         << std::endl;
      auto ov_indexer = [&](size_t iOb) -> int {
          return iOb * ov.nvars() + var_it_dist;
      };
      for (size_t sIndx = 0; sIndx < coords.n_obs; ++sIndx) {
        for (size_t rOb = 0; rOb < coords.record_counts[sIndx]; ++rOb) {
          const size_t obIdx = coords.record_starts[sIndx] + rOb;
          if (ov[ov_indexer(obIdx)] == missing_value) {
            data[obIdx] = typeToFill::value<T>();
          } else {
            data[obIdx] = ov[ov_indexer(obIdx)];
          }
        }
      }
  }
  void write_surf(NemoFeedbackWriter<T>& fdbk_writer,
      NemoFeedbackReduce& reducer, const std::string& nemo_name) {
      oops::Log::trace() << "NemoFeedback::variableData::write_surf "
                          << nemo_name << std::endl;
      reduced = reducer.reduce_data(data);
      fdbk_writer.write_variable_surf(nemo_name, reduced);
  }
  void write_profile(NemoFeedbackWriter<T>& fdbk_writer,
      NemoFeedbackReduce& reducer, const std::string& nemo_name) {
      oops::Log::trace() << "NemoFeedback::variableData::write_profile "
                         << nemo_name << std::endl;
      reducer.reduce_profile_data(data, reduced);
      const size_t n_reduced_profs = reducer.reduced_counts.size();
      if ((reducer.reduced_counts.at(n_reduced_profs-1)
           + reducer.reduced_starts.at(n_reduced_profs-1) > reduced.size()) ||
          (n_reduced_profs != reducer.n_obs_)) {
          std::ostringstream err_stream;
          err_stream << "NemoFeedback::variableData::write_profile extent "
                     << " data size " << data.size()
                     << " reduced size " << reduced.size()
                     << " n_reduced_profs " << n_reduced_profs
                     << " reducer.n_obs_ " << reducer.n_obs_
                     << " " << nemo_name << std::endl;
          err_stream << "NemoFeedback::variableData::write_profile "
                     << "index range out of bounds '" << nemo_name << "' "
                     << reducer.reduced_counts.at(n_reduced_profs-1) << " + "
                     << reducer.reduced_starts.at(n_reduced_profs-1) << " = "
                     << reducer.reduced_counts.at(n_reduced_profs-1) +
                      + reducer.reduced_starts.at(n_reduced_profs-1)
                     << " >= " << reduced.size();
          throw eckit::BadValue(err_stream.str(), Here());
      }
      fdbk_writer.write_variable_profile(nemo_name, reduced);
  }
};


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
    nameMap_(params.geoVaLsAliasFile.value()),
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
  const oops::Variables obsVarnames(varnames, channels);
  geovars_ = nameMap_.convertName(obsVarnames);
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

  if (obsdb_.comm().size() > 1) {
    std::stringstream ss;
    ss << (test_data_path.dirName() / test_data_path.baseName(false)).asString()
       << "_" << std::setw(5) << std::setfill('0') << obsdb_.comm().rank()
       << test_data_path.extension();
    test_data_path = eckit::PathName(ss.str());
  }

  // Generate lists of the variable names to setup the file
  NameData name_data;
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
    name_data.variable_names.push_back(nemoVariableParams.nemoName.value());
    name_data.legacy_ops_qc_conventions.push_back(obsdb_.has("QCFlags",
          nemoVariableParams.name.value()));
    name_data.long_names.push_back(nemoVariableParams.longName.value());
    name_data.unit_names.push_back(nemoVariableParams.units.value());
    extra_vars.push_back(nemoVariableParams.extravar.value().value_or(false));
    auto additionalVariablesParams = nemoVariableParams.variables.value();
    for (const NemoFeedbackAddVariableParameters& addVariableParams :
        additionalVariablesParams) {
      auto add_suffix = addVariableParams.feedbackSuffix.value();
      if (std::find(name_data.additional_names.begin(),
                    name_data.additional_names.end(), add_suffix)
          == name_data.additional_names.end()) {
        name_data.additional_names.push_back(add_suffix);
      }
    }
  }

  if (is_profile && is_altimeter)
    throw eckit::BadValue(std::string("NemoFeedback::postFilter cannot write")
        + " profile and altimeter data to the same file", Here());

  // obsdb_.nlocs is the number of point-observation locations on this
  // processor.
  size_t n_locs = obsdb_.nlocs();

  //  Calculate n_obs, n_levels, starts, counts
  //  Fill out lats, lons, depths, julian_days

  CoordData unreduced_coords;
  unreduced_coords.n_obs = n_locs;
  unreduced_coords.n_levels = 1;

  // Handle the where option.
  std::vector<bool> to_write = ufo::processWhere(parameters_.where, data_,
                                   ufo::WhereOperator::AND);
  if (to_write.size() != n_locs) to_write.assign(n_locs, true);

  groupCoordsByRecord(to_write, unreduced_coords, is_profile);

  // sync n_levels and juld_reference across files if multi-processing
  mpi_sync_coordinates(unreduced_coords, obsdb_.comm());

  // Calculate total number of obs to actually write.
  // This is already handled in the construction of record_counts/sizes above
  // for profiles, but for surface fields n_obs == n_locs, and we can look at
  // to_write
  const size_t n_surf_obs_to_write = is_profile ? unreduced_coords.n_obs
      : std::count(to_write.begin(), to_write.end(), true);

  // Set up the station type and identifier variables. The way to do this
  // depends on the data type. Also modify which obs to write if this is
  // altimeter data.
  std::vector<std::string> station_types(unreduced_coords.n_obs, "    ");
  std::vector<std::string> station_ids(unreduced_coords.n_obs, "        ");

  if (is_altimeter) {
    setupAltimeterIds(unreduced_coords.n_obs, station_ids, station_types, to_write);
  } else {
    setupIds(unreduced_coords.n_obs, unreduced_coords.record_starts, unreduced_coords.record_counts,
        station_ids, station_types);
  }

  NemoFeedbackReduce reducer(unreduced_coords.n_obs, n_surf_obs_to_write, to_write,
                             unreduced_coords.record_starts, unreduced_coords.record_counts);

  CoordData coords(std::move(unreduced_coords));
  if (is_profile) {
    coords.record_starts = reducer.reduced_starts;
    coords.record_counts = reducer.reduced_counts;
    std::vector<double> reduced_depths;
    reducer.reduce_profile_data(coords.depths, reduced_depths);
    coords.depths = reduced_depths;
  } else {
    coords.depths = reducer.reduce_data(coords.depths);
    coords.lats = reducer.reduce_data(coords.lats);
    coords.lons = reducer.reduce_data(coords.lons);
    coords.julian_days = reducer.reduce_data(coords.julian_days);
    station_ids = reducer.reduce_data(station_ids);
    station_types = reducer.reduce_data(station_types);
    coords.n_obs = n_surf_obs_to_write;
  }

  if (reducer.reduced_counts.size() != coords.n_obs) {
      std::ostringstream err_stream;
      err_stream << "nemo_feedback::NemoFeedback:: "
                 << "coords.n_obs does not match with dimension of reduced counts "
                 << " reduced_counts.size() " << reducer.reduced_counts.size()
                 << " coords.record_counts.size() " << coords.record_counts.size()
                 << " coords.n_obs " << coords.n_obs;
      throw eckit::BadValue(err_stream.str(), Here());
  }
  if (coords.record_counts.at(coords.n_obs-1)
      + coords.record_starts.at(coords.n_obs-1) > coords.depths.size()) {
      std::ostringstream err_stream;
      err_stream << "nemo_feedback::NemoFeedback:: "
                 << "index range out of bounds coords.n_obs "
                 << coords.n_obs << " record_starts.size() "
                 << coords.record_starts.size() << " counts: "
                 << coords.record_counts.at(coords.n_obs-1) << " + "
                 << coords.record_starts.at(coords.n_obs-1) << " = "
                 << coords.record_counts.at(coords.n_obs-1) +
                  + coords.record_starts.at(coords.n_obs-1)
                 << " >= " << coords.depths.size();
      throw eckit::BadValue(err_stream.str(), Here());
  }

  OutputDtype dtype = parameters_.type.value().value_or(OutputDtype::Double);
  if (dtype == OutputDtype::Float) {
    NemoFeedbackWriter<float> fdbk_writer(test_data_path,
                                          coords,
                                          name_data,
                                          extra_vars,
                                          station_types,
                                          station_ids);
    write_all_data<float> (fdbk_writer,
                            reducer,
                            coords,
                            ov,
                            is_profile);
  } else {
    NemoFeedbackWriter<double> fdbk_writer(test_data_path,
                                           coords,
                                           name_data,
                                           extra_vars,
                                           station_types,
                                           station_ids);
    write_all_data<double> (fdbk_writer,
                            reducer,
                            coords,
                            ov,
                            is_profile);
  }
}

template <typename T>
void NemoFeedback::write_all_data(NemoFeedbackWriter<T>& fdbk_writer,
                                  NemoFeedbackReduce& reducer,
                                  const CoordData& coords,
                                  const ioda::ObsVector &ov,
                                  const bool is_profile) const {
  // Write the data
  VariableData<T> variable_data{};

  size_t n_locs = obsdb_.nlocs();
  std::vector<int> variable_qcFlags(n_locs);
  std::vector<int> variable_qc(n_locs);
  std::vector<ufo::DiagnosticFlag> do_not_assimilate;
  for (const NemoFeedbackVariableParameters& nemoVariableParams :
        parameters_.variables.value()) {
    auto nemo_name = nemoVariableParams.nemoName.value();
    auto ufo_name = nemoVariableParams.name.value();
    auto obs_group = nemoVariableParams.iodaObsGroup.value()
        .value_or("ObsValue");
    variable_data.read_obs(obsdb_, obs_group, ufo_name);
    if (variable_data.data.size() != n_locs)
          throw eckit::BadValue("NemoFeedback::postFilter : no data ");
    auto extra_var = nemoVariableParams.extravar.value().value_or(false);
    if (!is_profile && extra_var) {
      variable_data.write_surf(fdbk_writer, reducer, nemo_name);
      // If this is an extra variable we do not want to write any of the
      // other variables with _OBS, _QC etc. added on to the name.
      continue;
    }

    if (is_profile) {
      variable_data.write_profile(fdbk_writer, reducer, nemo_name + "_OBS");
    } else {
      variable_data.write_surf(fdbk_writer, reducer, nemo_name + "_OBS");
    }

    // Write Met Office QC flag data for this variable if they exist.
    bool hasMetOfficeQC = false;
    if (obsdb_.has("QCFlags", ufo_name)) {
      obsdb_.get_db("QCFlags", ufo_name, variable_qcFlags);
      hasMetOfficeQC = true;
    } else {
      const size_t iv = flags_->varnames().find(ufo_name);
      variable_qcFlags.assign((*flags_)[iv].begin(), (*flags_)[iv].end());
    }

    // To the first qc flag index location
    std::vector<int> reduced_qcFlags;
    if (is_profile) {
      reducer.reduce_profile_data(variable_qcFlags, reduced_qcFlags);
      fdbk_writer.write_variable_level_qc(
          nemo_name + "_LEVEL_QC_FLAGS", reduced_qcFlags, 0);
    } else {
      std::vector<int> reduced_qcFlags = reducer.reduce_data(
          variable_qcFlags);
      fdbk_writer.write_variable_surf_qc(
          nemo_name + "_QC_FLAGS", reduced_qcFlags, 0);
    }

    // Whole Observation report QC flags
    if (is_profile) {
      std::vector<int> reduced_qc(coords.n_obs, 0);
      for (size_t iProf=0; iProf < coords.n_obs; ++iProf) {
        size_t badObs = 0;
        for (size_t iOb=0; iOb < coords.record_counts[iProf]; ++iOb) {
          if (hasMetOfficeQC &&
              (reduced_qcFlags[coords.record_starts[iProf] + iOb]
               & ufo::MetOfficeQCFlags::WholeObReport::FinalRejectReport)) {
            ++badObs;
          }
          if (!hasMetOfficeQC &&
              reduced_qcFlags[coords.record_starts[iProf] + iOb]) {
            ++badObs;
          }
        }
        reduced_qc[iProf] = (badObs == coords.record_counts[iProf] ? 4 : 1);
      }
      if (coords.record_counts.size() != coords.n_obs) {
        throw eckit::BadValue("NemoFeedback::postFilter : sizes don't match "
            + std::to_string(coords.record_counts.size()) + " != "
            + std::to_string(coords.n_obs), Here());
      }
      fdbk_writer.write_variable_surf_qc("OBSERVATION_QC", reduced_qc);
    } else {
      for (int i=0; i < variable_qc.size(); ++i) {
        if (hasMetOfficeQC &&
            (variable_qcFlags[i] &
             ufo::MetOfficeQCFlags::WholeObReport::FinalRejectReport)) {
          variable_qc[i] = 4;
        } else if (!hasMetOfficeQC && variable_qcFlags[i]) {
          variable_qc[i] = 4;
        } else {variable_qc[i] = 1;}
      }
      std::vector<int> reduced_qc = reducer.reduce_data(variable_qc);
      fdbk_writer.write_variable_surf_qc("OBSERVATION_QC", reduced_qc);
    }

    // Overall quality control flags
    if (obsdb_.has("DiagnosticFlags/FinalReject", ufo_name)) {
      std::vector<ufo::DiagnosticFlag> final_qc;
      obsdb_.get_db("DiagnosticFlags/FinalReject", ufo_name, final_qc);

      std::vector<int> variable_level_qc(n_locs);
      for (size_t iLoc = 0; iLoc < n_locs; ++iLoc) {
        if (final_qc[iLoc]) {
          variable_level_qc[iLoc] = 4;
        } else {variable_level_qc[iLoc] = 1;}
      }
      for (size_t iStart = 0; iStart < reducer.unreduced_starts.size();
           ++iStart) {
        size_t jLoc = reducer.unreduced_starts[iStart];
        if (final_qc[jLoc]) {
          variable_qc[iStart] = 4;
        } else {variable_qc[iStart] = 1;}
      }
      // Add do not assimilate flag if required.
      if (obsdb_.has("DiagnosticFlags/DoNotAssimilate", ufo_name)) {
        obsdb_.get_db("DiagnosticFlags/DoNotAssimilate", ufo_name,
            do_not_assimilate);
        for (size_t iLoc = 0; iLoc < n_locs; ++iLoc) {
          if (do_not_assimilate[iLoc]) {
            variable_level_qc[iLoc] += 128;
          }
        }
        for (size_t iStart = 0; iStart < reducer.unreduced_starts.size();
             ++iStart) {
          size_t jLoc = reducer.unreduced_starts[iStart];
          if (do_not_assimilate[jLoc]) {
            variable_qc[iStart] += 128;
          }
        }
      }
      std::vector<int> reduced_qc = reducer.reduce_data(variable_qc);
      fdbk_writer.write_variable_surf_qc(nemo_name + "_QC", reduced_qc);
      if (is_profile) {
        reducer.reduce_profile_data(variable_level_qc, reduced_qc);
        fdbk_writer.write_variable_level_qc(nemo_name + "_LEVEL_QC",
            reduced_qc);
      } else {
        reduced_qc = reducer.reduce_data(variable_level_qc);
        fdbk_writer.write_variable_surf_qc(nemo_name + "_LEVEL_QC", reduced_qc);
      }
    }

    // Write additional variables for this variable
    auto additionalVariablesParams = nemoVariableParams.variables.value();
    for (const NemoFeedbackAddVariableParameters& addParams :
        additionalVariablesParams) {
      auto add_name = nemo_name + "_" + addParams.feedbackSuffix.value();
      std::string ioda_group = addParams.iodaGroup.value();
      if (ioda_group == "HofX") {
        oops::Log::trace() << "NemoFeedback::HofX " << ov << std::endl;
        if (!ov.has(ufo_name)) {
          throw eckit::BadValue(
              "NemoFeedback::postFilter missing HofX variable: "
              + ioda_group + "/" + ufo_name, Here());
        }
        variable_data.read_hofx(ov, coords, ufo_name);
        if (is_profile) {
          variable_data.write_profile(fdbk_writer, reducer, add_name);
        } else {
          variable_data.write_surf(fdbk_writer, reducer, add_name);
        }
      } else if (!obsdb_.has(ioda_group, ufo_name)) {
        oops::Log::trace() << obsdb_ << std::endl;
        throw eckit::BadValue("NemoFeedback::postFilter missing variable: "
            + ioda_group + "/" + ufo_name, Here());
      } else {
        variable_data.read_obs(obsdb_, ioda_group, ufo_name);
        if (is_profile) {
          variable_data.write_profile(fdbk_writer, reducer, add_name);
        } else {
          variable_data.write_surf(fdbk_writer, reducer, add_name);
        }
      }
    }
  }
}

void NemoFeedback::groupCoordsByRecord(const std::vector<bool>& to_write,
                                       CoordData& coords,
                                       bool is_profile) const {
    coords.n_locs = obsdb_.nlocs();
    std::vector<util::DateTime> datetimes(coords.n_locs);
    obsdb_.get_db("MetaData", "dateTime", datetimes);
    coords.lats.resize(coords.n_locs);
    obsdb_.get_db("MetaData", "latitude", coords.lats);
    coords.lons.resize(coords.n_locs);
    obsdb_.get_db("MetaData", "longitude", coords.lons);

    if (coords.n_locs == 0) {
      coords.record_counts = std::vector<size_t>{0};
      coords.record_starts = std::vector<size_t>{0};
      coords.juld_reference = parameters_.refDate.value().value_or(
          util::DateTime{"1950-01-01T00:00:00Z"});
      return;
    }

    // Handle the reference date option.
    coords.juld_reference = parameters_.refDate.value().value_or(datetimes[0]);

    if (is_profile) {
      // if to_write all true set prune_profiles false
      bool all_obs_valid =
           std::all_of(to_write.begin(), to_write.end(),
               [](bool v) { return v; });
      bool prune_profiles = !all_obs_valid;

      std::vector<size_t> recnums = obsdb_.recidx_all_recnums();
      coords.n_levels = 0;
      std::vector<double> record_lats;
      std::vector<double> record_lons;
      std::vector<util::DateTime> record_dts;
      // guess average number of depths per profile to lessen reallocation
      size_t reclen_guess = 10;
      coords.record_counts.reserve(coords.n_locs/reclen_guess);
      coords.record_starts.reserve(coords.n_locs/reclen_guess);
      record_lats.reserve(coords.n_locs/reclen_guess);
      record_lons.reserve(coords.n_locs/reclen_guess);
      record_dts.reserve(coords.n_locs/reclen_guess);
      for (size_t iprof : recnums) {
        const std::vector<size_t> & obs_indices = obsdb_.recidx_vector(iprof);
        size_t n_levels_prof = 0;
        size_t reclen = 0;
        for (size_t jobs : obs_indices) {
          ++reclen;
          if (prune_profiles) {
            if (to_write[jobs]) { ++n_levels_prof; }
          } else {
            ++n_levels_prof;
          }
        }
        if (n_levels_prof != 0) {
          // we cannot guarantee profiles are stored in increasing depth order
          // (i.e "sort order: ascending") so we calculate the likely start
          // with a hack:
          size_t start = *std::min_element(obs_indices.begin(),
                                           obs_indices.end());
          coords.record_starts.push_back(start);
          coords.record_counts.push_back(reclen);
          record_lats.push_back(coords.lats[start]);
          record_lons.push_back(coords.lons[start]);
          record_dts.push_back(datetimes[start]);
        }
        if (n_levels_prof > coords.n_levels)
          coords.n_levels = n_levels_prof;
      }

      coords.n_obs = coords.record_counts.size();
      if (!prune_profiles && (recnums.size() != coords.n_obs)) {
        throw eckit::BadValue(std::string("NemoFeedback::groupCoordsByRecord ")
            + "recnums.size() != record_counts.size()",
            Here());
      }

      // n_levels is equal to largest number of levels across all the
      // profile data, provided we keep all the data in every profile.
      size_t n_levels_check = *std::max_element(coords.record_counts.begin(),
          coords.record_counts.end());
      if (!prune_profiles && (n_levels_check != coords.n_levels)) {
          throw eckit::BadValue(
              std::string("NemoFeedback::groupCoordsByRecord ")
              + "n_levels_check != coords.n_levels "
              + std::to_string(n_levels_check) + " "
              + std::to_string(coords.n_levels),
              Here());
      }

      coords.lats.resize(coords.n_obs);
      coords.lats.assign(record_lats.begin(), record_lats.end());
      coords.lons.resize(coords.n_obs);
      coords.lons.assign(record_lons.begin(), record_lons.end());
      datetimes.resize(coords.n_obs);
      datetimes.assign(record_dts.begin(), record_dts.end());
      coords.depths.resize(coords.n_locs);
      obsdb_.get_db(parameters_.depthGroup.value().value_or("MetaData"),
          parameters_.depthVariable.value().value_or("depthBelowWaterSurface"),
          coords.depths);
    } else {
      coords.n_obs = coords.n_locs;
      coords.depths.resize(coords.n_locs);
      std::fill(coords.depths.begin(), coords.depths.end(), 0);
      coords.record_counts.assign(coords.n_locs, 1);
      coords.record_starts.resize(coords.n_locs);
      for (int iLoc = 0; iLoc < coords.n_locs; ++iLoc)
        coords.record_starts[iLoc] = iLoc;
    }

    coords.julian_days.resize(coords.n_obs);
    for (int i=0; i < coords.n_obs; ++i) {
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
    // Define the station type variable. This may not always be defined in
    // obsdb_.
    if (obsdb_.has("MetaData", "fdbk_station_type")) {
      std::vector<int> station_types_int(obsdb_.nlocs());
      obsdb_.get_db("MetaData", "fdbk_station_type", station_types_int);
      char buffer[5];
      for (int iOb = 0; iOb < n_obs; ++iOb) {
        int jLoc = record_starts[iOb];
        snprintf(buffer, sizeof(buffer), "%4d", station_types_int[jLoc]);
        station_types[iOb] = buffer;
      }
    }

    // Define the station identifier variable. These may not always be defined
    // in odbsb_ and may be contained in either the station_id or
    // buoy_identifier variables, or both.
    if (obsdb_.has("MetaData", "stationIdentification")) {
      std::vector<std::string> station_ids_tmp(obsdb_.nlocs());
      obsdb_.get_db("MetaData", "stationIdentification", station_ids_tmp);
      std::string station_id_missing_value = util::missingValue(station_id_missing_value);
      for (int iOb = 0; iOb < n_obs; ++iOb) {
        int jLoc = record_starts[iOb];
        if ((station_ids_tmp[jLoc] != "") &&
            (station_ids_tmp[jLoc] != station_id_missing_value)) {
            station_ids_tmp[jLoc].resize(8, ' ');
            station_ids[iOb] = station_ids_tmp[jLoc].substr(0, 8);
        }
      }
    }
    if (obsdb_.has("MetaData", "buoyIdentifier")) {
      std::vector<int> buoy_ids(obsdb_.nlocs());
      obsdb_.get_db("MetaData", "buoyIdentifier", buoy_ids);
      char buffer[9];
      int buoy_id_missing_value = util::missingValue(buoy_id_missing_value);
      for (int iOb = 0; iOb < n_obs; ++iOb) {
        int jLoc = record_starts[iOb];
        if ((buoy_ids[jLoc] != buoy_id_missing_value) &&
            (station_ids[iOb] == std::string(8, ' '))){
          snprintf(buffer, sizeof(buffer), "%8d", buoy_ids[jLoc]);
          station_ids[iOb] = buffer;
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
    obsdb_.get_db("MetaData", "satelliteIdentifier", satellite_ids);

    // Get the most recent versions of the data.
    // Read version information.
    std::vector<int> version(n_obs);
    obsdb_.get_db("MetaData", "instrumentIdentifier", version);
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
  void NemoFeedback::mpi_sync_coordinates(CoordData& coords,
                                          const eckit::mpi::Comm& comm) {
  if (comm.size() > 0) {
    std::vector<size_t> all_nlevs(comm.size(), 0);
    comm.allGather(coords.n_levels, all_nlevs.begin(),
        all_nlevs.end());
    size_t max_nlev = *std::max_element(all_nlevs.begin(), all_nlevs.end());
    if (coords.n_levels < max_nlev) { coords.n_levels = max_nlev; }

    std::vector<double> juld_ref;
    coords.juld_reference.serialize(juld_ref);
    comm.broadcast(juld_ref.begin(), juld_ref.end(), 0);
    size_t index_0 = 0;
    coords.juld_reference.deserialize(juld_ref, index_0);
  }
}

void NemoFeedback::print(std::ostream & os) const {
  os << "NemoFeedback: config = " << parameters_ << std::endl;
}

}  // namespace nemo_feedback
