/*
 * (C) Crown Copyright 2020, the Met Office. All rights reserved.
 *
 * Refer to COPYRIGHT.txt of this distribution for details.
 */

#include "nemo-feedback/NemoFeedbackDataCreator.h"

#include <string.h>

#include <tuple>
#include <vector>
#include <memory>
#include <sstream>
#include <utility>

#include "eckit/exception/Exceptions.h"

#include "ioda/ObsDataVector.h"
#include "ioda/ObsSpace.h"
#include "ioda/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/DateTime.h"
#include "oops/util/missingValues.h"

#include "nemo-feedback/feedback_io/Utils.h"

namespace nemo_feedback {

NemoFeedbackDataCreator::NemoFeedbackDataCreator(
    const ioda::ObsSpace& obsdb, const ioda::ObsVector& hofxObsVector) :
  obsdb_(obsdb), hofxObsVector_(hofxObsVector) {
  size_t nLocations = obsdb.nlocs();

  std::vector<size_t> profileIndices = obsdb.recidx_all_recnums();

  std::vector<size_t> indices;
  std::vector<size_t> starts;

  indices.reserve(nLocations);
  starts.reserve(profileIndices.size());

  size_t nLevels = 0;
  for (size_t iProfile : profileIndices) {
    const std::vector<size_t> & obs_indices = obsdb.recidx_vector(iProfile);
    starts.push_back(indices.size());
    indices.insert(indices.end(), obs_indices.begin(), obs_indices.end());
    if (obs_indices.size() > nLevels)
      nLevels = obs_indices.size();
  }
  size_t nObs = profileIndices.size();

  indexer_ = std::make_shared<feedback_io::DataIndexer>(
      nObs, nLevels, nLocations, std::move(starts), std::move(indices));
}

template<typename T>
feedback_io::Data<T> NemoFeedbackDataCreator::create_from_obsdb(
    const std::string& obsGroup, const std::string& ufoName,
    const T typeInstance) const {
  oops::Log::trace() << NemoFeedbackDataCreator::className()
                     << ":create_from_obsdb " << obsGroup
                     << "/" << ufoName << std::endl;
  std::vector<T> data;
  obsdb_.get_db(obsGroup, ufoName, data);
  if (data.size() != obsdb_.nlocs())
    throw eckit::BadValue(NemoFeedbackDataCreator::className()
        + ":create_from_obsdb no data ");
  const T missingValue = util::missingValue<T>();
  // Convert oops missing values to NEMO missing value
  for_each(data.begin(), data.end(),
      [missingValue](T& d){ if (d == missingValue) {
                                d = feedback_io::typeToFill::value<T>();} });

  ASSERT_MSG(indexer_->n_locations() <= obsdb_.nlocs(),
      NemoFeedbackDataCreator::className()
      + ":create_from_obsdb indexer.n_locations() > obsdb_.nlocs()");
  return feedback_io::Data<T>(indexer_, std::move(data));
}

/// \brief create station ID feedbackData from obsdb
feedback_io::Data<std::string> NemoFeedbackDataCreator::create_from_obsdb(const
    std::string& obsGroup, const std::string& ufoName,
    const std::string TypeInstance, size_t width, bool leftJustify) const {
  oops::Log::trace() << NemoFeedbackDataCreator::className()
                     << ":create_from_obsdb " << obsGroup << "/" << ufoName
                     << std::endl;

  ASSERT_MSG(!leftJustify,
      "cannot left justify string variables - spaces are meaningful");

  std::vector<std::string> sourceData(obsdb_.nlocs());
  obsdb_.get_db(obsGroup, ufoName, sourceData);
  const std::string missingValue = util::missingValue<std::string()>();
  const std::string missingValueOut =
    feedback_io::typeToFill::value<std::string>();
  ASSERT_MSG(missingValueOut.size() == width,
      "width and missing value sizes do not match");
  std::stringstream stream;
  stream << "%" << width << "d";
  std::string format(stream.str());
  std::vector<std::string> data(obsdb_.nlocs(), missingValueOut);
  for (size_t iOb = 0; iOb < obsdb_.nlocs(); ++iOb) {
    if (missingValue != sourceData[iOb]) {
      sourceData[iOb].resize(width, ' ');
      data[iOb] = sourceData[iOb].substr(0, width);
    }
  }
  return feedback_io::Data<std::string>(indexer_, std::move(data));
}

/// \brief create station type feedbackData from integers in obsdb
feedback_io::Data<std::string> NemoFeedbackDataCreator::create_from_obsdb(
    const std::string& obsGroup, const std::string& ufoName,
    const int32_t TypeInstance, size_t width, bool leftJustify) const {
  oops::Log::trace() << NemoFeedbackDataCreator::className()
                     << ":create_from_obsdb " << obsGroup << "/" << ufoName
                     << std::endl;
  std::vector<int32_t> sourceData(obsdb_.nlocs(), 0);
  obsdb_.get_db(obsGroup, ufoName, sourceData);
  const int32_t missingValue = util::missingValue<int32_t>();
  const std::string missingValueOut =
    feedback_io::typeToFill::value<std::string>();
  ASSERT_MSG(missingValueOut.size() >= width,
      "width and missing value sizes do not match");

  std::stringstream stream;
  stream << "%";
  if (leftJustify)
    stream << "-";
  stream << width << "d";
  std::string format(stream.str());
  std::vector<std::string> data(obsdb_.nlocs(), missingValueOut);
  auto buffer = std::make_unique<char[]>(width + 1);
  for (size_t iOb = 0; iOb < obsdb_.nlocs(); ++iOb) {
    if (missingValue == sourceData[iOb]) {
      data[iOb] = missingValueOut.substr(0, width);
    } else {
      snprintf(buffer.get(), width + 1, format.c_str(), sourceData[iOb]);
      data[iOb] = std::string(buffer.get());
    }
  }
  return feedback_io::Data<std::string>(indexer_, std::move(data));
}

feedback_io::Data<int32_t> NemoFeedbackDataCreator::create_from_obsdb(const
    std::string& obsGroup, const std::string& ufoName, const
    ufo::DiagnosticFlag TypeInstance, int32_t whenTrue,
    int32_t whenFalse) const {
  oops::Log::trace() << NemoFeedbackDataCreator::className()
                     << ":create_from_obsdb ufo::DiagnosticFlag "
                     << obsGroup << "/" << ufoName << std::endl;
  std::vector<int32_t> data;
  std::vector<ufo::DiagnosticFlag> flagData;
  obsdb_.get_db(obsGroup, ufoName, flagData);
  for (ufo::DiagnosticFlag flag : flagData) {
    data.push_back(flag ? whenTrue : whenFalse);
  }
  if (data.size() != obsdb_.nlocs())
    throw eckit::BadValue(NemoFeedbackDataCreator::className()
        + ":create_from_obsdb ufo::DiagnosticFlag no data.");

  return feedback_io::Data<int32_t>(indexer_, std::move(data));
}

template<typename T>
feedback_io::Data<T> NemoFeedbackDataCreator::create_from_hofx(
    const std::string& ufoName, const T typeInstance) const {
  oops::Log::trace() << NemoFeedbackDataCreator::className()
                     << "::create_from_hofx HofX/" << ufoName << std::endl;
  std::vector<std::string> varNames = hofxObsVector_.varnames().variables();
  auto var_it = std::find(varNames.begin(), varNames.end(),
      ufoName);
  std::size_t var_it_dist = static_cast<std::size_t>(
                            std::distance(varNames.begin(), var_it));
  oops::Log::debug() << NemoFeedbackDataCreator::className()
                     << "::create_from_hofx :iterator distance is "
                     << var_it_dist << " hofxObsVector_.nvars: "
                     << hofxObsVector_.nvars() << std::endl;
  // Convert oops missing values to NEMO missing value
  const double missingValue = util::missingValue<double>();
  auto hofxIndexer = [&](size_t iOb) -> int {
      return iOb * hofxObsVector_.nvars() + var_it_dist;
  };
  std::vector<T> data;
  data.reserve(indexer_->n_locations());
  ASSERT_MSG(indexer_->n_locations() == hofxObsVector_.nlocs(),
      NemoFeedbackDataCreator::className()
      + ":create_from_hofx indexer.n_locations() != hofxObsVector_.nlocs()");
  for (size_t iOb = 0; iOb < hofxObsVector_.nlocs(); ++iOb) {
    if (hofxObsVector_[hofxIndexer(iOb)] == missingValue) {
      data.push_back(feedback_io::typeToFill::value<T>());
    } else {
      data.push_back(hofxObsVector_[hofxIndexer(iOb)]);
    }
  }
  return feedback_io::Data<T>(indexer_, std::move(data));
}

std::tuple<std::string, feedback_io::Data<double>>
  NemoFeedbackDataCreator::create_datetimes(
    const std::string& obsGroup, const std::string& ufoName,
    const util::DateTime juldReferenceDT) const {
  oops::Log::trace() << NemoFeedbackDataCreator::className()
      << ":create_datetimes " << obsGroup << "/" << ufoName
      << " with datetime " << juldReferenceDT << std::endl;
  if (obsGroup == "hofx" || obsGroup == "HofX") {
    std::ostringstream err_stream;
    err_stream << NemoFeedbackDataCreator::className()
               << "::create_from_hofx cannot retrieve a datetime HofX value."
               << std::endl;
    throw eckit::BadValue(err_stream.str(), Here());
  }

  std::string juldReference;
  {
    int year, month, day, hour, minute, second;
    juldReferenceDT.toYYYYMMDDhhmmss(year, month, day, hour, minute, second);
    std::ostringstream ref_stream;
    ref_stream << std::setfill('0');
    ref_stream << std::setw(4) << year;
    ref_stream << std::setw(2) << month;
    ref_stream << std::setw(2) << day;
    ref_stream << std::setw(2) << hour;
    ref_stream << std::setw(2) << minute;
    ref_stream << std::setw(2) << second;
    juldReference = ref_stream.str();
  }

  if (obsdb_.nlocs() == 0) {
    return std::make_tuple<std::string, feedback_io::Data<double>>(
        std::move(juldReference),
        std::move(feedback_io::Data<double>()));
  }

  std::vector<util::DateTime> datetimes;
  obsdb_.get_db(obsGroup, ufoName, datetimes);

  std::vector<double> data;
  data.reserve(obsdb_.nlocs());

  ASSERT_MSG(indexer_->n_locations() <= obsdb_.nlocs(),
      NemoFeedbackDataCreator::className()
      + ":create_from_obsdb indexer.n_locations() > obsdb_.nlocs()");

  for (size_t iLoc = 0; iLoc < datetimes.size(); ++iLoc) {
    util::Duration duration = static_cast<util::DateTime>(datetimes[iLoc])
        - juldReferenceDT;
    data.push_back(duration.toSeconds() / 86400.0);
  }

  return std::make_tuple<std::string, feedback_io::Data<double>>(
      std::move(juldReference),
      std::move(feedback_io::Data<double>(indexer_, std::move(data))));
}

template<typename T>
feedback_io::Data<std::string> NemoFeedbackDataCreator::create(const
    std::string& obsGroup, const std::string& ufoName, const T typeInstance,
    size_t width, bool leftJustify) const {
  if (obsGroup == "hofx" || obsGroup == "HofX") {
    std::ostringstream err_stream;
    err_stream << NemoFeedbackDataCreator::className()
               << "::create_from_hofx cannot retrieve string data."
               << std::endl;
    throw eckit::BadValue(err_stream.str(), Here());
  }
  if (!obsdb_.has(obsGroup, ufoName)) {
    throw eckit::BadValue(NemoFeedbackDataCreator::className() +
        ": missing obs variable: " + obsGroup + "/" + ufoName, Here());
  }
  return create_from_obsdb(obsGroup, ufoName, typeInstance, width, leftJustify);
}

template feedback_io::Data<std::string> NemoFeedbackDataCreator::create(
    const std::string& obsGroup, const std::string& ufoName,
    const int32_t typeInstance, size_t width, bool leftJustify) const;
template feedback_io::Data<std::string> NemoFeedbackDataCreator::create(
    const std::string& obsGroup, const std::string& ufoName,
    const std::string typeInstance, size_t width, bool leftJustify) const;

std::tuple<feedback_io::Data<std::string>, feedback_io::Data<std::string>>
NemoFeedbackDataCreator::create_altimeter_IDs() const {
    // Station type and station identifier are both defined from the
    // satellite identifier metadata for Altimetry data
    std::vector<int32_t> satellite_ids(obsdb_.nlocs());
    std::vector<std::string> station_ids(obsdb_.nlocs());
    std::vector<std::string> station_types(obsdb_.nlocs());
    obsdb_.get_db("MetaData", "satelliteIdentifier", satellite_ids);

    char buffer9[9];
    char buffer5[5];
    for (size_t iLoc = 0; iLoc < obsdb_.nlocs(); ++iLoc) {
      snprintf(buffer9, sizeof(buffer9), "%04d    ", satellite_ids[iLoc]);
      station_ids[iLoc] = static_cast<std::string>(buffer9);
      snprintf(buffer5, sizeof(buffer5), "%4d", satellite_ids[iLoc]);
      station_types[iLoc] = static_cast<std::string>(buffer5);
    }

    return std::make_tuple<feedback_io::Data<std::string>,
                           feedback_io::Data<std::string>>(
        feedback_io::Data<std::string>(indexer_, std::move(station_ids)),
        feedback_io::Data<std::string>(indexer_, std::move(station_types)));
}

feedback_io::Data<int32_t> NemoFeedbackDataCreator::create(const std::string&
    obsGroup, const std::string& ufoName, const ufo::DiagnosticFlag
    typeInstance, int32_t whenTrue, int32_t whenFalse) const {
  if (obsGroup == "hofx" || obsGroup == "HofX") {
    std::ostringstream err_stream;
    err_stream << NemoFeedbackDataCreator::className()
               << "::create_from_hofx cannot retrieve a"
               << " ufo::DiagnosticFlag HofX value." << std::endl;
    throw eckit::BadValue(err_stream.str(), Here());
  }
  if (!obsdb_.has(obsGroup, ufoName)) {
    throw eckit::BadValue(NemoFeedbackDataCreator::className() +
        ": missing obs variable: " + obsGroup + "/" + ufoName, Here());
  }
  return create_from_obsdb(obsGroup, ufoName, typeInstance, whenTrue,
      whenFalse);
}

template<typename T>
feedback_io::Data<T> NemoFeedbackDataCreator::create(
    const std::string& obsGroup, const std::string& ufoName,
    const T typeInstance) const {
  if (obsGroup == "hofx" || obsGroup == "HofX") {
    if (!hofxObsVector_.has(ufoName)) {
      throw eckit::BadValue(NemoFeedbackDataCreator::className() +
          ": missing HofX variable: " + ufoName, Here());
    }
    return create_from_hofx(ufoName, typeInstance);
  }
  if (!obsdb_.has(obsGroup, ufoName)) {
    throw eckit::BadValue(NemoFeedbackDataCreator::className() +
        ": missing obs variable: " + obsGroup + "/" + ufoName, Here());
  }
  return create_from_obsdb(obsGroup, ufoName, typeInstance);
}

template feedback_io::Data<int32_t> NemoFeedbackDataCreator::create(
    const std::string& obsGroup, const std::string& ufoName,
    const int32_t typeInstance) const;
template feedback_io::Data<double> NemoFeedbackDataCreator::create(
    const std::string& obsGroup, const std::string& ufoName,
    const double typeInstance) const;
template feedback_io::Data<float> NemoFeedbackDataCreator::create(
    const std::string& obsGroup, const std::string& ufoName,
    const float typeInstance) const;

}  // namespace nemo_feedback
