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

#include "sgd.h"
#include <cmath>

namespace fdb {
namespace learning {

L2SSVM::L2SSVM(double l, double reg)
    : learning_rate_(l), regularization_constant_(reg) {
  iteration_ = 0;
  w_.clear();
}

void L2SSVM::UpdateLearningRate() {
  iteration_++;
  learning_rate_ /= iteration_ * .01 + 1;
}

std::tuple<double, double> L2SSVM::Update(const FeatureVec& true_fv,
                                          const FeatureVec& wmv_fv,
                                          const double cur_loss) {
  const double true_score = CalculateSparseScalarProduct(w_, true_fv);
  const double wmv_score = CalculateSparseScalarProduct(w_, wmv_fv);
  const double slack_value = wmv_score + cur_loss - true_score;
  if (slack_value > 0.0) {
    FeatureVec gradient;
    AddSparseMultiple(-1.0, true_fv, &gradient);
    AddSparseMultiple(1.0, wmv_fv, &gradient);
    // Update model according to gradient information.
    for (auto gradient_it = gradient.cbegin(); gradient_it != gradient.cend();
         ++gradient_it) {
      const int feature_id = gradient_it->first;
      const double feature_value = gradient_it->second;
      // Add feature contribution to the score of this example.
      double gradient_value = learning_rate_ * feature_value;
      auto w_it = w_.find(feature_id);
      if (w_it != w_.end()) {
        w_it->second -=
            regularization_constant_ * w_it->second + gradient_value;
      } else {
        w_.insert(make_pair(feature_id, gradient_value));
      }
    }
  }
  return std::make_tuple(max(0.0, slack_value), wmv_score);
}

// Implementation of the L1-regularized SSVM.
L1SSVM::L1SSVM(double l, double reg, int num_ex)
    : learning_rate_(l), regularization_constant_(reg), num_examples_(num_ex) {
  iteration_ = 0;
  alpha_ = 0.85;
  accumulated_l1_penalty_ = 0.0;
  w_.clear();
  q_.clear();
}

double L1SSVM::LearningRateAt(int k) {
  double cur_learning_rate =
      pow(alpha_, -static_cast<double>(k) / num_examples_);
  cur_learning_rate *= learning_rate_;
  return cur_learning_rate;
}

std::tuple<double, double> L1SSVM::Update(const FeatureVec& true_fv,
                                          const FeatureVec& wmv_fv,
                                          const double cur_loss) {
  const double true_score = CalculateSparseScalarProduct(w_, true_fv);
  const double wmv_score = CalculateSparseScalarProduct(w_, wmv_fv);
  double slack_value = wmv_score + cur_loss - true_score;
  if (slack_value > 0.0) {
    FeatureVec gradient;
    AddSparseMultiple(-1.0, true_fv, &gradient);
    AddSparseMultiple(1.0, wmv_fv, &gradient);
    double cur_learning_rate = LearningRateAt(iteration_);
    accumulated_l1_penalty_ +=
        cur_learning_rate * regularization_constant_ / num_examples_;
    // Update model according to gradient information.
    for (auto gradient_it = gradient.cbegin(); gradient_it != gradient.cend();
         ++gradient_it) {
      const int feature_id = gradient_it->first;
      const double feature_value = gradient_it->second;
      // Add feature contribution to the score of this example.
      auto w_it = w_.find(feature_id);
      auto q_it = q_.find(feature_id);
      double qi = 0.0;
      if (q_it != q_.end()) {
        qi = q_it->second;
      }

      double wi = cur_learning_rate * -feature_value;
      if (w_it != w_.end()) {
        wi = w_it->second + cur_learning_rate * -feature_value;
      }

      double z = wi;
      if (wi > 0.0) wi = max(0.0, wi - (accumulated_l1_penalty_ + qi));
      if (wi < 0.0) wi = min(0.0, wi + (accumulated_l1_penalty_ - qi));

      auto insert_ret = w_.insert(make_pair(feature_id, wi));
      if (!insert_ret.second) insert_ret.first->second = wi;

      insert_ret = q_.insert(make_pair(feature_id, wi - z));
      if (!insert_ret.second) insert_ret.first->second += wi - z;
    }
    iteration_++;
  }
  return std::make_tuple(std::max(0.0, slack_value), wmv_score);
}

}  // namespace learning
}  // namespace fdb
