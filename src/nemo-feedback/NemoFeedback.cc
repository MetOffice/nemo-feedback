/*
 * (C) Crown Copyright 2020, the Met Office. All rights reserved.
 *
 * Refer to COPYRIGHT.txt of this distribution for details. 
 */

#include "nemo-feedback/NemoFeedback.h"

#include <utility>
#include <vector>
#include <bitset>

#include "eckit/config/Configuration.h"
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
#include "ufo/utils/metoffice/MetOfficeQCFlags.h"


namespace nemo_feedback {

NemoFeedback::NemoFeedback(ioda::ObsSpace & obsdb, const Parameters_ & params,
                           std::shared_ptr<ioda::ObsDataVector<int> > flags,
                           std::shared_ptr<ioda::ObsDataVector<float> > obsErrors)
  : obsdb_(obsdb), geovars_(), flags_(std::move(flags)), obsErrors_(std::move(obsErrors)),
    parameters_(params)
{
  oops::Log::trace() << "NemoFeedback constructor starting" << std::endl;

  const std::vector<int> channels{};
  std::vector<std::string> varnames;
  for (const eckit::LocalConfiguration& nemo_variable_conf : parameters_.Variables.value()) {
    const std::string varname = nemo_variable_conf.getString("name");
    if(std::find(varnames.begin(), varnames.end(), varname) == varnames.end())
      varnames.push_back(varname); 
  }
  geovars_ = oops::Variables(varnames, channels);
  //extradiagvars_ = oops::Variables(varnames, channels); 

  eckit::LocalConfiguration conf(parameters_.toConfiguration());
  // Validity time is set to the midpoint of the assimilation window
  const util::DateTime validityTime =
      obsdb.windowStart() + (obsdb.windowEnd() - obsdb.windowStart()) / 2;
  conf.set("validity_time", validityTime.toString());

  MPI_Comm mpiComm = MPI_COMM_WORLD;
  if (auto parallelComm = dynamic_cast<const eckit::mpi::Parallel*>(&obsdb.comm())) {
    mpiComm = parallelComm->MPIComm();
  }

}

NemoFeedback::~NemoFeedback() {
  oops::Log::trace() << "NemoFeedback destructor " << std::endl;

}

void NemoFeedback::priorFilter(const ufo::GeoVaLs & gv) {
  oops::Log::trace() << "NemoFeedback priorFilter" << std::endl;

}

//void NemoFeedback::postFilter(const ioda::ObsVector &ov,
//                              const ufo::ObsDiagnostics &dv) {
void NemoFeedback::postFilter(const ioda::ObsVector &ov,
                  const ioda::ObsVector &bv,
                  const ufo::ObsDiagnostics &dv) {
  oops::Log::trace() << "NemoFeedback postFilter" << std::endl;

  eckit::PathName test_data_path(parameters_.Filename);
  // ov.n_obs is the number of non-missing non-qced obs?
  // obsdb_.nlocs is the number of original locations I think.
  size_t n_obs = obsdb_.nlocs();
  size_t n_levels = 1;
  std::vector<double> lats(n_obs,0);
  obsdb_.get_db("MetaData", "latitude", lats);
  std::vector<double> lons(n_obs,0);
  obsdb_.get_db("MetaData", "longitude", lons);
  std::vector<double> depths(n_obs*n_levels,0);

  std::vector<std::string> datetimes(n_obs);
  obsdb_.get_db("MetaData", "datetime", datetimes);

  util::DateTime juld_reference(datetimes[0]);

  std::vector<double> julian_days(n_obs,0);
  for (int i=0; i<n_obs; ++i) {
    julian_days[i] = (static_cast<util::DateTime>(datetimes[i]) - juld_reference).toSeconds() /
      86400.0;
  }

  // Generate lists of the variable names to setup the file
  std::vector<std::string> additional_names;
  std::vector<std::string> variable_names;
  for (const eckit::LocalConfiguration& nemo_variable_conf : parameters_.Variables.value()) {
    variable_names.push_back(nemo_variable_conf.getString("nemo name"));
    std::vector<eckit::LocalConfiguration> additional_variables;
    nemo_variable_conf.get("additional variables", additional_variables); 
    for (const eckit::LocalConfiguration& var_conf : additional_variables) {
      auto add_suffix = var_conf.getString("feedback suffix");
      if (std::find(additional_names.begin(), additional_names.end(), add_suffix) == additional_names.end()) {
        additional_names.push_back(add_suffix);
      }
    }
  }

  NemoFeedbackWriter fdbk_writer( test_data_path, lons, lats, depths,
      julian_days, variable_names, additional_names, n_levels, juld_reference);

  // Write the data
  std::vector<double> variable_data;
  std::vector<int> variable_qcFlags;
  std::vector<int> variable_qc;
  for (const eckit::LocalConfiguration& nemo_variable_conf : parameters_.Variables.value()) {

    auto nemo_name = nemo_variable_conf.getString("nemo name");
    auto ufo_name = nemo_variable_conf.getString("name");
    obsdb_.get_db("ObsValue", ufo_name, variable_data);
    auto missing_value = util::missingValue(variable_data[0]);
    for(int i=0; i<ov.size(); ++i) {
      if (variable_data[i] == missing_value) variable_data[i] = NemoFeedbackWriter::double_fillvalue;
    }
    fdbk_writer.write_variable_surf(nemo_name + "_OBS", variable_data);

    // Write QC flag data for this variable to the first qc flag index location
    obsdb_.get_db("QCFlags", ufo_name, variable_qcFlags);
    fdbk_writer.write_variable_surf_qc(nemo_name + "_QC_FLAGS", variable_qcFlags, 0);

    // Convert Met Office QC flags to Ocean Quality Control flags
    variable_qc.resize(variable_qcFlags.size());
    for (int i=0; i<variable_qc.size(); ++i) {
      if(variable_qcFlags[i] & ufo::MetOfficeQCFlags::Elem::BackRejectFlag) {
        variable_qc[i]=4;
      } else {variable_qc[i]=0;}
    }
    fdbk_writer.write_variable_surf_qc(nemo_name + "_QC", variable_qc);

    // Whole Observation report QC flags
    for (int i=0; i<variable_qc.size(); ++i) {
      if(variable_qcFlags[i] & ufo::MetOfficeQCFlags::WholeObReport::FinalRejectReport) {
        variable_qc[i]=4;
      } else {variable_qc[i]=0;}
    }
    fdbk_writer.write_variable_surf_qc("OBSERVATION_QC", variable_qc);

    // Write additional variables for this variable
    std::vector<eckit::LocalConfiguration> additional_variables;
    nemo_variable_conf.get("additional variables", additional_variables); 
    for (const eckit::LocalConfiguration& var_conf : additional_variables) {
      auto add_name = nemo_name + "_" + var_conf.getString("feedback suffix");
      auto ioda_group = var_conf.getString("ioda group");
      if (ioda_group == "HofX") {
        oops::Log::trace() << ov << std::endl;
        if (ov.has(ufo_name)) {
          std::vector<std::string> ov_varnames = ov.varnames().variables();
          auto var_it = std::find(ov_varnames.begin(), ov_varnames.end(), ufo_name);
          oops::Log::trace() << "ov_varnames = " << ov_varnames << std::endl;
          oops::Log::trace() << "iterator distance is "
            << static_cast<std::size_t>(std::distance(ov_varnames.begin(), var_it))
            << std::endl;
          auto missing_value_add = util::missingValue(ov[0]);
          oops::Log::trace() << "Missing value: " << missing_value_add << std::endl;
          for(int i=0; i<ov.size(); ++i) {
            if (ov[i] == missing_value_add) {
              variable_data[i] = NemoFeedbackWriter::double_fillvalue;
            } else {
              variable_data[i] = ov[i];
            }
          }
          fdbk_writer.write_variable_surf(add_name, variable_data);
        }
      } else if (!obsdb_.has(ioda_group, ufo_name)) {
        oops::Log::trace() << obsdb_ << std::endl;
        throw eckit::BadValue("missing variable: " + ioda_group + "/" + ufo_name, Here());
      } else {
        obsdb_.get_db(ioda_group, ufo_name, variable_data);
        fdbk_writer.write_variable_surf(add_name, variable_data);
      }
    }
  }

}

void NemoFeedback::print(std::ostream & os) const {
  os << "NemoFeedback: config = " << parameters_ << std::endl;
}

}  // namespace nemo_feedback
