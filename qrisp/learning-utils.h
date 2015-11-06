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

#ifndef FDB_LEARNING_SGD_UTILS_H_
#define FDB_LEARNING_SGD_UTILS_H_

#include "utils.h"

#include <unordered_map>

namespace fdb {
namespace learning {

using namespace std;
using qrisp::FeatureVec;

// This projection to the l1-ball was implemented using the description in
// "Efficient Projections onto the l1-Ball for Learning in High-Dimensions".
//
void perform_l1_projection(FeatureVec* cur_w);

// This projection to the l2-ball was implemented using the description in
// "Pegasos: Primal Estimated sub-GrAdient SOlver for SVM"
void perform_l2_projection(double lambda, FeatureVec* cur_w);

// This function calculates the following formula:
// w[i] += scalar * v[i]
inline void AddSparseMultiple(double scalar, const FeatureVec& v,
                              FeatureVec* w) {
  for (auto it = v.cbegin(); it != v.cend(); ++it) {
    // Add feature contribution to the score of this example.
    const double& feature_value = it->second;
    auto insert_ret = w->insert(make_pair(it->first, scalar * feature_value));
    if (!insert_ret.second) {
      insert_ret.first->second += scalar * feature_value;
    }
  }
}

inline void CalculateSparseScalarMultiplication(const double& scalar,
                                                FeatureVec* v) {
  for (auto it = v->begin(); it != v->end(); ++it) {
    it->second *= scalar;
  }
}

inline double SparseScalarProductAux(const FeatureVec& smaller,
                                     const FeatureVec& larger) {
  double scalar_prod = 0.0;
  for (auto it = smaller.cbegin(); it != smaller.cend(); ++it) {
    // Add feature contribution to the score of this example.
    FeatureVec::const_iterator larger_entry = larger.find(it->first);
    if (larger_entry != larger.end()) {
      scalar_prod += larger_entry->second * it->second;
    }
  }
  return scalar_prod;
}

inline double CalculateSparseScalarProduct(const FeatureVec& w,
                                           const FeatureVec& v) {
  if (w.size() <= v.size()) {
    return SparseScalarProductAux(w, v);
  }
  return SparseScalarProductAux(v, w);
}

}  // namespace learning
}  // namespace fdb

#endif  // FDB_LEARNING_SGD_UTILS_H_
