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

#ifndef EXPERIMENTAL_USERS_FDB_LEARNING_SGD_SGD_H_
#define EXPERIMENTAL_USERS_FDB_LEARNING_SGD_SGD_H_

#include "learning-utils.h"

#include <algorithm>
#include <iterator>

namespace fdb {
namespace learning {

// This implementation is based on...
//
class L2SSVM {
 public:
  L2SSVM(const double l, const double reg);

  std::tuple<double, double> Update(const FeatureVec& true_fv,
                                    const FeatureVec& wmv_fv,
                                    const double cur_loss);

  void UpdateLearningRate();

  FeatureVec w_;
  double learning_rate_;

 private:
  double regularization_constant_;
  int iteration_;
};

// This implementation is based on the paper
// "Stochastic Gradient Descent Training for L1-regularized Log-linear Models
// with Cumulative Penalty"
// Y. Tsuruoka,  J. Tsujii, S. Ananiadou
// Proceedings of the 47th Annual Meeting of the ACL and the 4th IJCNLP of the
// AFNLP, pages 477–485,
class L1SSVM {
 public:
  L1SSVM(const double l, const double reg, const int num_ex);

  std::tuple<double, double> Update(const FeatureVec& true_fv,
                                    const FeatureVec& wmv_fv,
                                    const double cur_loss);

  FeatureVec w_;

 private:
  double LearningRateAt(int k);

  double learning_rate_;
  double regularization_constant_;
  double alpha_;
  double accumulated_l1_penalty_;
  int iteration_;
  int num_examples_;
  FeatureVec q_;
};

}  // namespace learning
}  // namespace fdb
#endif  // EXPERIMENTAL_USERS_FDB_LEARNING_SGD_SGD_H_
