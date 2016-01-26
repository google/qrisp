// Copyright 2014 Google Inc.All Rights Reserved.
//
//     Licensed under the Apache License,
//     Version 2.0(the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>

#include <google/protobuf/stubs/status.h>
#include <google/protobuf/stubs/status_macros.h>
#include <google/protobuf/stubs/statusor.h>

#include "dataset-utils.h"
#include "model.h"
#include "performance.h"
#include "qrisp-lib.h"
#include "qrisp/proto/config.pb.h"
#include "recurrences-nbest.h"
#include "recurrences.h"
#include "sgd.h"
#include "structure.h"
#include "third_party/googleflags/include/gflags/gflags.h"
#include "utils.h"

DEFINE_string(configuration_file, "",
              "Location of the asciipb representation of the "
              "ConfigMessage message.");

using namespace std;

using ::google::protobuf::util::Status;
using ::google::protobuf::util::error::INVALID_ARGUMENT;
using ::google::protobuf::util::error::OUT_OF_RANGE;

namespace {

// TODO(fdb): Some preliminary checks. Add more and use proto inspection instead
// of nasty "if blah return blub".
Status ValidateConfigMessage(const qrisp::ConfigMessage& config) {
  if (config.learning_rate() <= 0.0 || config.learning_rate() > 1.0) {
    return Status(OUT_OF_RANGE, " learning_rate. Valid range: (0.0, 1.0].");
  }
  if (config.regularization_coeff() <= 0.0 ||
      config.regularization_coeff() > 10000.0) {
    return Status(OUT_OF_RANGE,
                  " regularization_coeff. Valid range: (0.0, 10000.0].");
  }
  if (config.max_iter() > 1000) {
    return Status(OUT_OF_RANGE, " max_iter. Valid range: (-inf, 1000].");
  }
  return Status();
}

}  // namespace

int main(int argc, char** argv) {
  google::InitGoogleLogging(argv[0]);
  google::ParseCommandLineFlags(&argc, &argv, true);

  if (FLAGS_configuration_file.empty()) {
    LOG(ERROR) << "Found empty path for configuration file.";
    return EXIT_FAILURE;
  }

  qrisp::ConfigMessage config;
  LoadProtoFromFile(FLAGS_configuration_file, &config);

  LOG(INFO) << "Validating ConfigMessage..";
  const auto status = ValidateConfigMessage(config);
  if (!status.ok()) {
    LOG(ERROR) << status;
    return EXIT_FAILURE;
  }

  // Load structures.
  qrisp::Dataset data;
  qrisp::LoadDataset(config.input_data(), &data);
  if (data.size() == 0) {
    LOG(ERROR) << "Empty dataset.";
    return EXIT_FAILURE;
  }
  LOG(INFO) << "Loaded data. Size is: " << data.size();

  // Subselect training data.
  qrisp::Dataset working_set;
  if (!config.split_fn().empty()) {
    vector<string> split;
    qrisp::LoadSplit(config.split_fn(), &split);
    qrisp::SubselectData(split, data, &working_set);
    LOG(INFO) << "Filtered working_set. New size is: " << working_set.size();
  } else {
    std::swap(data, working_set);
  }

  // Load holdout set.
  qrisp::Dataset holdout_set;
  if (!config.holdout_fn().empty()) {
    vector<string> split;
    qrisp::LoadSplit(config.holdout_fn(), &split);
    qrisp::SubselectData(split, data, &holdout_set);
    LOG(INFO) << "Filtered holdout set. New size is: " << holdout_set.size();
  }

  // Load existing model if given.
  qrisp::QRSPModel model;
  if (!config.pretrained_model_fn().empty()) {
    qrisp::LoadProtoFromFile(config.pretrained_model_fn(), &model);
  }

  LOG(INFO) << "max_iter is: " << config.max_iter() << " starting in mode: ";
  if (config.max_iter() > 0) {
    LOG(INFO) << "-- Training";
    qrisp::StartTraining(working_set, holdout_set, model, config);
  } else {
    LOG(INFO) << "-- Prediction";
    qrisp::StartPrediction(working_set, model, config);
  }
  return EXIT_SUCCESS;
}
