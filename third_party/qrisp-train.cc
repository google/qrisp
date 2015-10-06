// This program is part of the QRSP project. It is used to train a model.

#include <algorithm>
#include <cmath>
#include <string>
#include <sstream>

#include "dataset-utils.h"
#include "model.h"
#include "performance.h"
#include "proto/config.pb.h"
#include "recurrences-nbest.h"
#include "recurrences.h"
#include "rna-structure.h"
#include "sgd.h"
#include "utils.h"
#include <gflags/gflags.h>

DECLARE_bool(enable_quality_features);

DEFINE_bool(plot_features, false,
            "Output all feature vectors of training set.");

DEFINE_bool(visualize_training, false,
            "Show what is going on during training.");

DEFINE_bool(resume_training, true, "Uses stored training model, if present.");

DEFINE_double(epsilon, 0.0001, "The precision used during training.");

DEFINE_double(learning_rate, 1.0, "");

DEFINE_double(regularization_coeff, 1.0, "");

DEFINE_string(holdout_fn, "", "Location of the split file.");

DEFINE_string(resume_model_fn, "", "Location of the split file.");

DEFINE_string(input_sstable, "", "");

DEFINE_string(pretrained_model_fn, "",
              "A file containing a pre-trained model.");

DEFINE_string(result_dir, "", "A file containing a pre-trained model.");

DEFINE_string(split_fn, "", "Location of the split file.");

DEFINE_int32(max_iter, 1, "Max. number of passes through the whole dataset.");

DEFINE_string(nbest_values, "1",
              "List of values for nbest-centroid "
              "calculation. Example: If the list is '1,3,7' the centroids from "
              "the 1, 3 and 7 highest-scoring structures are used as the "
              "predicted structures. So basically we report performance @1, "
              "@3 and @7.");

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

void StartPrediction(const Dataset& dataset, const QRSPModel& model,
                     TrainingConfiguration* cfg) {
  LOG(INFO) << "Starting prediction...";
  vector<int> cluster_sizes;
  SplitStringToInts(FLAGS_nbest_values, &cluster_sizes);
  std::sort(cluster_sizes.begin(), cluster_sizes.end());
  cluster_sizes.erase(std::unique(cluster_sizes.begin(), cluster_sizes.end()),
                      cluster_sizes.end());
  const int nbest = cluster_sizes.back();
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
                   const QRSPModel& model, TrainingConfiguration* cfg) {
  CHECK(training_set.size() > 0) << "Empty training set.";
  LOG(INFO) << "Starting training...";
  // Extract different cluster sizes and make values unique.
  vector<int> cluster_sizes(cfg->nbest_values().begin(),
                            cfg->nbest_values().end());
  std::sort(cluster_sizes.begin(), cluster_sizes.end());
  cluster_sizes.erase(std::unique(cluster_sizes.begin(), cluster_sizes.end()),
                      cluster_sizes.end());
  CHECK(cluster_sizes.size() > 0);
  CHECK(cluster_sizes.back() <= nbest::max_nbest);
  LOG(INFO) << "Largest cluster size: " << cluster_sizes.back();

  const idx_t num_examples = training_set.size();
  const double epsilon = FLAGS_epsilon;
  // Variables shared by all modes.
  int max_slack = 0;
  for (const auto& elem : training_set) {
    max_slack += elem.second.GetSize();
  }
  LOG(INFO) << "Max. theoretic slack: " << max_slack;
  // Initialize stochastic gradient descent solver.
  // L1SSVM ssvm(cfg->learning_rate(), cfg->regularization_coeff(),
  // num_examples);
  L2SSVM ssvm(cfg->learning_rate(), cfg->regularization_coeff());
  VLOG(1) << "Learning rate / regularization: " << cfg->learning_rate() << " "
          << cfg->regularization_coeff();
  FeatureVec cur_w;
  if (model.feature_size() > 0) {
    InitializeFeatureVecFromModel(model, &cur_w);
  }
  int i = 1;
  if (FLAGS_resume_training) {
    QRSPModel model;
    LOG(INFO) << "Using existing model: " << FLAGS_resume_model_fn;
    LOG(INFO) << "Resuming training at iteration: " << i;
    LoadProtoFromFile<QRSPModel>(FLAGS_resume_model_fn, &model);
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
  for (; i <= cfg->max_iter(); i++) {
    double sum_of_slacks = 0.0;
    // Inner loop over examples.
    std::shuffle(ids.begin(), ids.end(), std::default_random_engine(s));
    for (const auto& key : ids) {
      const Structure& true_structure = training_set.find(key)->second;
      if (FLAGS_enable_quality_features && !true_structure.HasQuality()) {
        LOG(ERROR) << "Quality features enabled, but structure has no quality "
                      "information.";
      }
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
    oss << FLAGS_result_dir << "/model.iter_" << i;
    SaveProtoToFile<QRSPModel>(oss.str(), resulting_model);
    const double obj_value =
        0.5 * sqrt(CalculateSparseScalarProduct(ssvm.w_, ssvm.w_)) +
        sum_of_slacks * cfg->regularization_coeff();
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
  SaveProtoToFile<QRSPModel>(FLAGS_result_dir + "/final_model.txt",
                             resulting_model);
}

}  // namespace qrisp

namespace {

void SetConfigurationFromFlags(qrisp::TrainingConfiguration* cfg) {
  cfg->set_enable_quality_features(FLAGS_enable_quality_features);
  cfg->set_plot_features(FLAGS_plot_features);
  cfg->set_visualize_training(FLAGS_visualize_training);
  //
  cfg->set_epsilon(FLAGS_epsilon);
  cfg->set_learning_rate(FLAGS_learning_rate);
  cfg->set_regularization_coeff(FLAGS_regularization_coeff);
  //
  cfg->set_pretrained_model_fn(FLAGS_pretrained_model_fn);
  cfg->set_split_fn(FLAGS_split_fn);
  cfg->set_holdout_fn(FLAGS_holdout_fn);
  cfg->set_input_sstable(FLAGS_input_sstable);
  cfg->set_result_dir(FLAGS_result_dir);
  //
  cfg->set_max_iter(FLAGS_max_iter);

  vector<int> values;
  qrisp::SplitStringToInts(FLAGS_nbest_values, &values);
  cfg->mutable_nbest_values()->Resize(values.size(), 0);
  std::copy(values.begin(), values.end(), cfg->mutable_nbest_values()->begin());
  CHECK(cfg->nbest_values_size() > 0)
      << "Error: At least one value required for field: nbest_values";
}

}  // namespace

int main(int argc, char** argv) {
  google::InitGoogleLogging(argv[0]);

  // Check if command line settings make sense at all.
  CHECK(FLAGS_learning_rate > 0.0 && FLAGS_learning_rate <= 100.0)
      << "Invalid learning rate...";
  CHECK(FLAGS_regularization_coeff > 0.0 &&
        FLAGS_regularization_coeff <= 10000.0)
      << "Invalid regularization factor...";
  CHECK(!FLAGS_input_sstable.empty())
      << "You have to supply a file with training data.";
  CHECK(FLAGS_max_iter >= -1 && FLAGS_max_iter <= 1000)
      << "Invalid range for max. number of iterations (-1 to 1.000).";
  // Load structures.
  qrisp::Dataset data;
  qrisp::LoadDataset(FLAGS_input_sstable, &data);
  CHECK(data.size() > 0) << "Empty dataset.";
  LOG(INFO) << "Loaded data. Size is: " << data.size();
  // Subselect training data.
  qrisp::Dataset working_set;
  if (!FLAGS_split_fn.empty()) {
    vector<string> split;
    qrisp::LoadSplit(FLAGS_split_fn, &split);
    qrisp::SubselectData(split, data, &working_set);
    LOG(INFO) << "Filtered working_set. New size is: " << working_set.size();
  } else {
    std::swap(data, working_set);
  }
  qrisp::Dataset holdout_set;
  if (!FLAGS_holdout_fn.empty()) {
    vector<string> split;
    qrisp::LoadSplit(FLAGS_holdout_fn, &split);
    qrisp::SubselectData(split, data, &holdout_set);
    LOG(INFO) << "Filtered holdout set. New size is: " << holdout_set.size();
  }
  // Save current setup.
  qrisp::TrainingConfiguration cfg;
  SetConfigurationFromFlags(&cfg);
  qrisp::SaveProtoToFile<qrisp::TrainingConfiguration>(
      FLAGS_result_dir + "/config.txt", cfg);
  // Load existing model if given.
  qrisp::QRSPModel model;
  if (!FLAGS_pretrained_model_fn.empty()) {
    qrisp::LoadProtoFromFile<qrisp::QRSPModel>(FLAGS_pretrained_model_fn,
                                               &model);
  }
  if (FLAGS_max_iter > 0) {
    qrisp::StartTraining(working_set, holdout_set, model, &cfg);
  } else {
    qrisp::StartPrediction(working_set, model, &cfg);
  }
  return EXIT_SUCCESS;
}
