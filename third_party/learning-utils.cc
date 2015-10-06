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

#include "learning-utils.h"

#include <algorithm>
#include <vector>

namespace fdb {
namespace learning {

using namespace std;

void perform_l1_projection(FeatureVec* cur_w) {
  double z = 1.0;
  vector<double> mu(cur_w->size(), 0.0);
  size_t jdx = 0;
  for (auto w_it = cur_w->begin(); w_it != cur_w->end(); ++w_it) {
    mu[jdx] = w_it->second;
    jdx++;
  }
  std::sort(mu.begin(), mu.end());
  std::reverse(mu.begin(), mu.end());
  size_t rho = 0;
  double cum_sum = 0.0;
  double max_value = 0.0;
  double value = 0.0;
  for (size_t idx = 0; idx < mu.size(); idx++) {
    cum_sum += mu[idx] - z;
    value = mu[idx] - cum_sum / static_cast<double>(idx + 1);
    if (value > max_value) {
      max_value = value;
      rho = idx;
    }
  }
  double theta = 0.0;
  for (size_t idx = 0; idx <= rho; idx++) theta += mu[idx] - z;
  theta /= static_cast<double>(rho + 1);

  for (auto w_it = cur_w->begin(); w_it != cur_w->end(); ++w_it)
    w_it->second = std::max(w_it->second - theta, 0.0);
}

void perform_l2_projection(double lambda, FeatureVec* cur_w) {
  double scalar_prod = CalculateSparseScalarProduct(*cur_w, *cur_w);
  double factor = min(1.0, (1.0 / sqrt(lambda)) / sqrt(scalar_prod));
  FeatureVec::iterator w_it;
  for (w_it = cur_w->begin(); w_it != cur_w->end(); ++w_it)
    w_it->second *= factor;
}

}  // namespace learning
}  // namespace fdb
