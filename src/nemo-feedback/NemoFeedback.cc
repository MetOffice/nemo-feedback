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
#include "nemo-feedback/NemoFeedbackDataCreator.h"
#include "nemo-feedback/feedback_io/Utils.h"
#include "nemo-feedback/feedback_io/Writer.h"
#include "nemo-feedback/feedback_io/Data.h"
#include "nemo-feedback/feedback_io/DataIndexer.h"
#include "ufo/GeoVaLs.h"
#include "ufo/ObsDiagnostics.h"
#include "ufo/filters/DiagnosticFlag.h"
#include "ufo/filters/processWhere.h"
#include "ufo/filters/ObsAccessor.h"
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

  // Generate lists of the variable names to setup the file
  bool isProfile = false;
  isAltimeter_ = false;
  for (const NemoFeedbackVariableParameters& nemoVariableParams :
        parameters_.variables.value()) {
    if (nemoVariableParams.nemoName.value() == "SLA") {
      isAltimeter_ = true;
    }
    if (nemoVariableParams.nemoName.value() == "POTM" ||
        nemoVariableParams.nemoName.value() == "PSAL") {
      isProfile = true;
    }
    nameData_.variable_names.push_back(nemoVariableParams.nemoName.value());
    nameData_.legacy_ops_qc_conventions.push_back(obsdb_.has("QCFlags",
          nemoVariableParams.name.value()));
    nameData_.long_names.push_back(nemoVariableParams.longName.value());
    nameData_.unit_names.push_back(nemoVariableParams.units.value());
    isExtraVariable_.push_back(nemoVariableParams.extravar.value()
        .value_or(false));
    auto additionalVariablesParams = nemoVariableParams.variables.value();
    for (const NemoFeedbackAddVariableParameters& addVariableParams :
        additionalVariablesParams) {
      auto add_suffix = addVariableParams.feedbackSuffix.value();
      if (std::find(nameData_.additional_names.begin(),
                    nameData_.additional_names.end(), add_suffix)
          == nameData_.additional_names.end()) {
        nameData_.additional_names.push_back(add_suffix);
      }
    }
  }

  if (isProfile && isAltimeter_)
    throw eckit::BadValue(std::string("NemoFeedback::postFilter cannot write")
        + " profile and altimeter data to the same file", Here());
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

  // Handle the where option.
  std::vector<bool> to_write = ufo::processWhere(parameters_.where, data_,
                                   ufo::WhereOperator::AND);
  if (to_write.size() != obsdb_.nlocs()) to_write.assign(obsdb_.nlocs(), true);

  // exclude all but the latest altimetry observations
  if (isAltimeter_) {
    updateAltimeterSelection(to_write);
  }

  NemoFeedbackDataCreator creator(obsdb_, ov, to_write);

  feedback_io::MetaData metaData(setupMetaData(creator));

  OutputDtype dtype = parameters_.type.value().value_or(OutputDtype::Double);
  if (dtype == OutputDtype::Float) {
    feedback_io::Writer<float> fdbk_writer(test_data_path,
                                           metaData,
                                           nameData_,
                                           isExtraVariable_);
    write_all_data<float> (fdbk_writer, creator);
  } else {
    feedback_io::Writer<double> fdbk_writer(test_data_path,
                                            metaData,
                                            nameData_,
                                            isExtraVariable_);
    write_all_data<double> (fdbk_writer, creator);
  }
  oops::Log::trace() << "NemoFeedback postFilter done" << std::endl;
}

template <typename T>
void NemoFeedback::write_all_data(feedback_io::Writer<T>& fdbk_writer,
                                  const NemoFeedbackDataCreator& creator) const
{
  std::vector<ufo::DiagnosticFlag> do_not_assimilate;
  for (const NemoFeedbackVariableParameters& nemoVariableParams :
        parameters_.variables.value()) {
    auto nemo_name = nemoVariableParams.nemoName.value();
    auto ufo_name = nemoVariableParams.name.value();
    auto obs_group = nemoVariableParams.iodaObsGroup.value()
        .value_or("ObsValue");

    feedback_io::Data<T> variableData(creator.create(obs_group, ufo_name,
          T(0)));

    auto extra_var = nemoVariableParams.extravar.value().value_or(false);
    if (extra_var) {
      fdbk_writer.write_variable_profile(nemo_name, variableData);
      // If this is an extra variable we do not want to write any of the
      // other variables with _OBS, _QC etc. added on to the name.
      continue;
    }

    fdbk_writer.write_variable_profile(nemo_name + "_OBS", variableData);

    // Write Met Office QC flag data for this variable if they exist.
    feedback_io::Data<int32_t> variableQCFlagsData;
    bool hasMetOfficeQC = false;
    if (obsdb_.has("QCFlags", ufo_name)) {
      variableQCFlagsData = creator.create("QCFlags", ufo_name, int32_t(0));
      hasMetOfficeQC = true;
    } else {
      const size_t iv = flags_->varnames().find(ufo_name);
      std::vector<int32_t> variable_qcFlags;
      variable_qcFlags.assign((*flags_)[iv].begin(), (*flags_)[iv].end());
      variableQCFlagsData = feedback_io::Data<int32_t>(creator.indexer(),
          variable_qcFlags);
    }

    // If profile data, write to level QC Flags
    if (variableData.n_levels() == 1) {
      fdbk_writer.write_variable_surf_qc(
          nemo_name + "_QC_FLAGS", variableQCFlagsData, 0);
    } else {
      fdbk_writer.write_variable_level_qc(
          nemo_name + "_LEVEL_QC_FLAGS", variableQCFlagsData, 0);
    }

    // Whole Observation report QC flags
    feedback_io::Data<int32_t> variableQCData(variableData.indexer(),
        std::vector<int32_t>(variableData.n_locations(), 1));
    {
      for (size_t iProf = 0; iProf < variableData.n_obs(); ++iProf) {
        size_t badObs = 0;
        for (size_t iLevel = 0; iLevel < variableData.length(iProf); ++iLevel) {
          if (hasMetOfficeQC &&
              (variableQCFlagsData(iProf, iLevel)
               & ufo::MetOfficeQCFlags::WholeObReport::FinalRejectReport)) {
            ++badObs;
          }
          if (!hasMetOfficeQC &&
              variableQCFlagsData(iProf, iLevel)) {
            ++badObs;
          }
        }
        variableQCData[iProf] = (badObs == variableData.length(iProf) ? 4 : 1);
      }
      fdbk_writer.write_variable_surf_qc("OBSERVATION_QC", variableQCData);
    }

    // Overall quality control flags
    feedback_io::Data<int32_t> variableFinalQCData;

    if (obsdb_.has("DiagnosticFlags/FinalReject", ufo_name)) {
      std::vector<ufo::DiagnosticFlag> final_qc;
      variableFinalQCData = creator.create("DiagnosticFlags/FinalReject",
          ufo_name, ufo::DiagnosticFlag(0), 4, 1);

      // Add do not assimilate flag if required.
      if (obsdb_.has("DiagnosticFlags/DoNotAssimilate", ufo_name)) {
        feedback_io::Data<int32_t> variableDoNotAssimilateData(
            creator.create("DiagnosticFlags/DoNotAssimilate", ufo_name,
              ufo::DiagnosticFlag(0), 128, 0));
        for (size_t iProfile = 0;
            iProfile < variableDoNotAssimilateData.n_obs(); ++iProfile) {
          for (size_t iLevel = 0;
              iLevel < variableDoNotAssimilateData.length(iProfile);
              ++iLevel) {
            variableFinalQCData(iProfile, iLevel) +=
              variableDoNotAssimilateData(iProfile, iLevel);
          }
        }
      }

      fdbk_writer.write_variable_surf_qc(nemo_name + "_QC",
          variableFinalQCData);
      fdbk_writer.write_variable_level_qc(nemo_name + "_LEVEL_QC",
          variableFinalQCData);
    }

    // Write additional variables for this variable
    auto additionalVariablesParams = nemoVariableParams.variables.value();
    for (const NemoFeedbackAddVariableParameters& addParams :
        additionalVariablesParams) {
      auto add_name = nemo_name + "_" + addParams.feedbackSuffix.value();
      std::string ioda_group = addParams.iodaGroup.value();
      feedback_io::Data<T> variableAdditionalData(creator.create(ioda_group,
            ufo_name, T(0)));
      fdbk_writer.write_variable_profile(add_name, variableAdditionalData);
    }
  }
}

feedback_io::MetaData NemoFeedback::setupMetaData(
    const NemoFeedbackDataCreator& creator) const {

  feedback_io::Data<double> lats(creator.create("MetaData", "latitude",
        static_cast<double>(0)));
  feedback_io::Data<double> lons(creator.create("MetaData", "longitude",
        static_cast<double>(0)));
  const std::string depthGroup = parameters_.depthGroup.value()
    .value_or("MetaData");
  const std::string depthVariable = parameters_.depthVariable.value()
    .value_or("depthBelowWaterSurface");
  feedback_io::Data<double> depths;
  if (obsdb_.has(depthGroup, depthVariable)) {
    depths = feedback_io::Data<double>(creator.create(depthGroup, depthVariable,
          static_cast<double>(0)));
  } else {
    depths = feedback_io::Data<double>(creator.indexer(),
        std::vector<double>(obsdb_.nlocs(), 0));
  }

  auto [juldReferenceGlobal, nLevelsGlobal] = mpiSync(lats.n_levels());  // NOLINT(*)

  auto [juldReference, julianDays] = creator.create_datetimes("MetaData",  // NOLINT(*)
      "dateTime", juldReferenceGlobal);

  feedback_io::Data<std::string> stationIDs, stationTypes;
  if (isAltimeter_) {
    std::tie(stationIDs, stationTypes) = creator.create_altimeter_IDs();
  } else {
    std::tie(stationIDs, stationTypes) = setupIDs(creator);
  }

  return feedback_io::MetaData(lats,
                  lons,
                  julianDays,
                  depths,
                  stationTypes,
                  stationIDs,
                  nLevelsGlobal,
                  juldReference);
}

std::tuple<feedback_io::Data<std::string>, feedback_io::Data<std::string>>
  NemoFeedback::setupIDs(const NemoFeedbackDataCreator& creator) const {
  oops::Log::trace() << "NemoFeedback::setupIDs: starting" << std::endl;
  constexpr size_t stationIDWidth = 8;
  feedback_io::Data<std::string> stationIDs;
  bool stationIdentificationAvailable = false;
  if (obsdb_.has("MetaData", "stationIdentification")) {
    stationIDs = feedback_io::Data<std::string>(creator.create("MetaData",
          "stationIdentification",
          std::string(""), stationIDWidth));
    stationIdentificationAvailable = true;
  }

  int32_t buoyIDMissingValueInt = util::missingValue(buoyIDMissingValueInt);
  constexpr size_t buoyIDWidth = 8;
  const std::string missingStringFeedback =
    feedback_io::typeToFill::value<std::string>();

  if (obsdb_.has("MetaData", "buoyIdentifier")) {
    feedback_io::Data<std::string> buoyIDs(creator.create("MetaData",
          "buoyIdentifier", int32_t(0), buoyIDWidth));
    for (size_t iOb = 0; iOb < stationIDs.n_obs(); ++iOb) {
      if (stationIdentificationAvailable) {
        if (buoyIDs[iOb] != missingStringFeedback &&
            stationIDs[iOb] == missingStringFeedback) {
          stationIDs[iOb] = buoyIDs[iOb];
        }
      } else {
          stationIDs[iOb] = buoyIDs[iOb];
      }
    }
    stationIdentificationAvailable = true;
  }

  if (!stationIdentificationAvailable) {
    std::vector<std::string> blankStationIDData(obsdb_.nlocs(),
        std::string(8, ' '));
    stationIDs = feedback_io::Data<std::string>(creator.indexer(),
        blankStationIDData);
  }

  feedback_io::Data<std::string> stationTypes;
  constexpr size_t stationTypeWidth = 8;
  if (obsdb_.has("MetaData", "fdbk_station_type")) {
    stationTypes = feedback_io::Data<std::string>(creator.create("MetaData",
          "fdbk_station_type", int32_t(0), stationTypeWidth));
  } else {
    std::vector<std::string> blankStationTypeData(obsdb_.nlocs(),
        std::string(4, ' '));
    stationTypes = feedback_io::Data<std::string>(creator.indexer(),
        blankStationTypeData);
  }

  return std::make_tuple<
         feedback_io::Data<std::string>,
         feedback_io::Data<std::string>> (
            std::move(stationIDs), std::move(stationTypes));
}

// Filter to retrieve the most recent version of the data
void NemoFeedback::updateAltimeterSelection(std::vector<bool>& to_write) const {
  size_t n_obs = obsdb_.nlocs();
  std::vector<util::DateTime> datetimes(n_obs);
  obsdb_.get_db("MetaData", "dateTime", datetimes);

  // Station type and station identifier are both defined from the
  // satellite_identifier metadata, so we use this to identify the
  // data source.
  std::vector<int> satellite_ids(n_obs);
  obsdb_.get_db("MetaData", "satelliteIdentifier", satellite_ids);

  std::vector<int> version(n_obs);
  obsdb_.get_db("MetaData", "instrumentIdentifier", version);

  // Get integer values for the day of each data point.
  std::vector<int> ymd(n_obs);
  int ymd_tmp;
  int hms_tmp;
  for (size_t i = 0; i < n_obs; ++i) {
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
      for (size_t i=0; i < n_obs; ++i) {
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
        for (size_t i=0; i < n_obs; ++i) {
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

std::tuple<util::DateTime, size_t> NemoFeedback::mpiSync(size_t nLevelsLocal)
  const {
  auto& comm = obsdb_.comm();
  util::DateTime juldReferenceLocal = parameters_.refDate.value()
    .value_or(util::DateTime{"1950-01-01T00:00:00Z"});
  std::vector<util::DateTime> datetimes;
  obsdb_.get_db("MetaData", "dateTime", datetimes);
  juldReferenceLocal = parameters_.refDate.value().value_or(datetimes[0]);
  oops::Log::trace() << "NemoFeedback::mpiSync " << juldReferenceLocal
                     << " and " << nLevelsLocal << " levels" << std::endl;

  util::DateTime juldReferenceGlobal(juldReferenceLocal);
  size_t nLevelsGlobal = nLevelsLocal;
  if (comm.size() == 0) {
    return std::make_tuple<util::DateTime, size_t>(
        std::move(juldReferenceGlobal),
        std::move(nLevelsGlobal));
  }

  std::vector<size_t> all_nlevs(comm.size(), 0);
  comm.allGather(nLevelsLocal, all_nlevs.begin(), all_nlevs.end());
  nLevelsGlobal = *std::max_element(all_nlevs.begin(), all_nlevs.end());

  if (!parameters_.refDate.value()) {
    std::vector<double> juldRef;
    juldReferenceLocal.serialize(juldRef);
    size_t index_0 = 0;
    comm.broadcast(juldRef.begin(), juldRef.end(), 0);
    juldReferenceGlobal.deserialize(juldRef, index_0);
  }

  return std::make_tuple<util::DateTime, size_t>(std::move(juldReferenceGlobal),
      std::move(nLevelsGlobal));
}

void NemoFeedback::print(std::ostream & os) const {
  os << "NemoFeedback: config = " << parameters_ << std::endl;
}
}  // namespace nemo_feedback
