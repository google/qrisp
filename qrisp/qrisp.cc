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

#include <gflags/gflags.h>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>

#include "dataset-utils.h"
#include "model.h"
#include "performance.h"
#include "qrisp-lib.h"
#include "qrisp/proto/config.pb.h"
#include "recurrences-nbest.h"
#include "recurrences.h"
#include "rna-structure.h"
#include "sgd.h"
#include "utils.h"

DEFINE_bool(resume_training, true, "Uses stored training model, if present.");
DEFINE_string(resume_model_fn, "", "Location of the split file.");

DEFINE_string(settings, "",
              "Location of the ascii representation of the "
              "Configuration message.");

using namespace std;

namespace {

// TODO(fdb): Some preliminary checks. Add more and use proto inspection instead
// of nasty "if blah return blub".
bool CheckConfiguration(const qrisp::Configuration& config) {
  if (!(config.learning_rate() > 0.0 && config.learning_rate() <= 100.0)) {
    return false;
  }
  if (!(config.regularization_coeff() > 0.0 &&
        config.regularization_coeff() <= 10000.0)) {
    return false;
  }
  if (!(config.max_iter() >= -1 && config.max_iter() <= 1000)) {
    return false;
  }
  return true;
}

}  // namespace

int main(int argc, char** argv) {
  google::InitGoogleLogging(argv[0]);

  if (FLAGS_settings.empty()) {
    LOG(ERROR) << "Found empty path for settings.";
    return EXIT_FAILURE;
  }

  qrisp::Configuration config;
  LoadProtoFromFile(FLAGS_settings, &config);

  // Load structures.
  qrisp::Dataset data;
  qrisp::LoadDataset(config.input_data(), &data);
  CHECK(data.size() > 0) << "Empty dataset.";
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
  if (config.pretrained_model_fn().empty()) {
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
