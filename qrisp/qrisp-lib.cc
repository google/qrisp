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
#include <random>
#include <sstream>
#include <string>

#include "model.h"
#include "performance.h"
#include "qrisp-lib.h"
#include "recurrences-nbest.h"
#include "recurrences.h"
#include "sgd.h"
#include "structure.h"
#include "utils.h"

using namespace std;

namespace qrisp {

using fdb::learning::CalculateSparseScalarProduct;
using fdb::learning::L1SSVM;
using fdb::learning::L2SSVM;

typedef const pair<string, Structure>& Instance;

Performance UpdateBestPerformance(const Performance& previous,
                                  const Performance& current) {
  if (F1Measure(previous) < F1Measure(current)) {
    return current;
  }
  return previous;
}

void CalculatePredictions(const FeatureVec& model, const Dataset& data,
                          const vector<int>& cluster_sizes,
                          Predictions* predictions) {
  // Elements of cluster_sizes are already sorted in increasing order.
  const int nbest = cluster_sizes.back();
  for (const auto& elem : data) {
    predictions->push_back(vector<vector<int>>());
    if (nbest > 1) {
      nbest::DecodePaths(nbest, elem.second, model, LOSS_DISABLED, 1.0,
                         &(predictions->back()));
    } else {
      vector<int> result;
      DecodeBestPath(elem.second, model, LOSS_DISABLED, 1.0, &result);
      predictions->back().push_back(result);
    }
  }
}

int ExtractClusterSizes(const ConfigMessage& config, vector<int>* sizes) {
  if (config.nbest_values_size() > 0) {
    sizes->assign(config.nbest_values().begin(), config.nbest_values().end());
    std::sort(sizes->begin(), sizes->end());
    std::sort(sizes->begin(), sizes->end());
    sizes->erase(std::unique(sizes->begin(), sizes->end()), sizes->end());
    return sizes->back();
  }
  sizes->assign(1, 1);
  return 1;
}

void StartPrediction(const Dataset& dataset, const QRSPModel& model,
                     const ConfigMessage& config) {
  vector<int> cluster_sizes;
  const int nbest = ExtractClusterSizes(config, &cluster_sizes);
  LOG(INFO) << "Largest cluster size: " << nbest;

  // Load model used in prediction.
  FeatureVec m;
  if (model.feature_size() > 0) {
    LOG(INFO) << "Loading model...";

    InitializeFeatureVecFromModel(model, &m);
  }
  // Elements of cluster_sizes are already sorted in increasing order.
  Predictions predictions;
  for (const auto& elem : dataset) {
    predictions.push_back(vector<vector<int>>());
    nbest::DecodePaths(nbest, elem.second, m, LOSS_DISABLED, 1.0,
                       &(predictions.back()));
  }
  // Perform predictions.
  PerformanceMetrics metrics;
  EstimatePerformanceOnHoldout(dataset, predictions, cluster_sizes, &metrics);
  for (const auto& elem : metrics) {
    LOG(INFO) << "Cluster size: " << elem.first;
    LOG(INFO) << "Sens. / spec.: " << elem.second.first << " "
              << elem.second.second;
  }
}

void StartTraining(const Dataset& training_set, const Dataset& holdout_set,
                   const QRSPModel& model, const ConfigMessage& config) {
  CHECK(training_set.size() > 0) << "Empty training set.";

  vector<int> cluster_sizes;
  const int nbest = ExtractClusterSizes(config, &cluster_sizes);
  LOG(INFO) << "Largest cluster size: " << nbest;

  const idx_t num_examples = training_set.size();
  const double epsilon = config.epsilon();
  // Variables shared by all modes.
  int max_slack = 0;
  for (const auto& elem : training_set) {
    max_slack += elem.second.Size();
  }
  LOG(INFO) << "Max. theoretic slack: " << max_slack;
  // Initialize stochastic gradient descent solver.
  // L1SSVM ssvm(config.learning_rate(), config.regularization_coeff(),
  // num_examples);
  L2SSVM ssvm(config.learning_rate(), config.regularization_coeff());
  VLOG(1) << "Learning rate / regularization: " << config.learning_rate() << " "
          << config.regularization_coeff();
  FeatureVec cur_w;
  if (model.feature_size() > 0) {
    InitializeFeatureVecFromModel(model, &cur_w);
  }
  int i = 1;
  if (config.resume_training()) {
    QRSPModel model;
    LOG(INFO) << "Using existing model: " << config.resume_model_fn();
    LOG(INFO) << "Resuming training at iteration: " << i;
    LoadProtoFromFile<QRSPModel>(config.resume_model_fn(), &model);
    InitializeFeatureVecFromModel(model, &cur_w);
  }
  ssvm.w_ = cur_w;
  LOG(INFO) << "Size of model vector: " << cur_w.size();
  // Generate a set of indices of the training data that can be shuffled at
  // every epoch. This is useful for stochastic gradient descent.
  vector<string> ids(num_examples, "");
  transform(training_set.begin(), training_set.end(), ids.begin(),
            [](Instance e) -> string { return e.first; });
  CHECK(ids.size() == training_set.size());
  const double s = 12.34;
  std::tuple<double, double> best_performance = std::make_pair(0.0, 0.0);
  LOG(INFO) << "Entering training loop..";
  // Outer loop over training epochs.
  for (; i <= config.max_iter(); i++) {
    double sum_of_slacks = 0.0;
    // Inner loop over examples.
    std::shuffle(ids.begin(), ids.end(), std::default_random_engine(s));
    for (const auto& key : ids) {
      const Structure& true_structure = training_set.find(key)->second;
      vector<idx_t> true_pairings;
      true_structure.GetPairings(&true_pairings);
      // Decode new structure based on model.
      vector<idx_t> decoded_pairings;
      const double dp_score = DecodeBestPath(
          true_structure, cur_w, LOSS_ENABLED, 1.0, &decoded_pairings);
      Structure wmv_structure(true_structure);
      wmv_structure.SetPairings(decoded_pairings);
      const double cur_loss = HammingLoss(true_structure, wmv_structure);
      FeatureVec wmv_fv;
      CalculateFeatures(wmv_structure, &wmv_fv);
      const double scalar_prod = CalculateSparseScalarProduct(cur_w, wmv_fv);
      FeatureVec true_fv;
      CalculateFeatures(true_structure, &true_fv);
      score_t slack;
      score_t wmv_score;
      std::tie(slack, wmv_score) = ssvm.Update(true_fv, wmv_fv, cur_loss);
      VLOG(1) << "Current loss / slack: " << cur_loss << " / " << slack;
      cur_w = ssvm.w_;
      sum_of_slacks += slack;
      if (fabs(dp_score - (wmv_score + cur_loss)) > epsilon ||
          fabs(dp_score - (scalar_prod + cur_loss)) > epsilon) {
        LOG(ERROR) << "Inconsistent scores: dp_score = " << dp_score
                   << " wmv_score + loss: " << wmv_score + cur_loss
                   << " scalar_prod + loss: " << scalar_prod + cur_loss;
      }
    }
    ssvm.UpdateLearningRate();
    QRSPModel resulting_model;
    // Now estimate performance on holdout set.
    if (!holdout_set.empty()) {
      Predictions predictions;
      CalculatePredictions(cur_w, holdout_set, cluster_sizes, &predictions);
      PerformanceMetrics metrics;
      EstimatePerformanceOnHoldout(holdout_set, predictions, cluster_sizes,
                                   &metrics);
      const auto it = metrics.find(1);
      if (it != metrics.end()) {
        auto current = std::make_tuple(it->second.first, it->second.second);
        best_performance = UpdateBestPerformance(best_performance, current);
        // DisplayHoldoutPerformance(i, best_performance, current);
      }
      for (const auto& m : metrics) {
        auto* perf = resulting_model.add_holdout_performance();
        perf->set_iteration(i);
        perf->set_cluster_size(m.first);
        perf->set_sensitivity(m.second.first);
        perf->set_specificity(m.second.second);
      }
    }
    InitializeModelFromFeatureVec(cur_w, &resulting_model);
    std::ostringstream oss;
    oss << config.result_dir() << "/model.iter_" << i;
    SaveProtoToFile<QRSPModel>(oss.str(), resulting_model);
    const double obj_value =
        0.5 * sqrt(CalculateSparseScalarProduct(ssvm.w_, ssvm.w_)) +
        sum_of_slacks * config.regularization_coeff();
    LOG(INFO) << "Iteration: " << i << ".\n"
              << " -- Objective: " << obj_value << "\n"
              << " -- Slacks (sum): " << sum_of_slacks << "\n"
              << " -- learning rate: " << ssvm.learning_rate_;
    if (!(sum_of_slacks > 0.0)) {
      LOG(INFO) << "Complete pass over training set w/o update. Stopping...";
      break;
    }
  }
  QRSPModel resulting_model;
  InitializeModelFromFeatureVec(cur_w, &resulting_model);
  SaveProtoToFile<QRSPModel>(config.result_dir() + "/final_model.txt",
                             resulting_model);
}

}  // namespace qrisp
