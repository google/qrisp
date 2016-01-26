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

#ifndef QRISP_MODEL_H_
#define QRISP_MODEL_H_

#include "structure.h"
#include "utils.h"

#include <vector>

namespace qrisp {

static const idx_t MAX_LOOP_SIZE = 30;

// Definition of function ptr for scoring functions.
typedef void (*ScoringFunction_t)(const Tuple& t, const Structure& rna,
                                  IndexVec* indices);

void UpdateFeatureWeights(const IndexVec& vec, FeatureVec* fv);

score_t GetSummedParameters(const IndexVec& vec, const FeatureVec& params);

// We can observe 10 different possible base pairings.
// AA AC AG AU
// CC CG CU
// GG GU
// UU
// 00 01 02 03 11 12 13 22 23 33
//  0  1  2  3  4  5  6  7  8  9
inline idx_t GetBasePairOffset(idx_t i, idx_t j) {
  if (i > 3 || j > 3) {
    exit(EXIT_FAILURE);
  }
  if (i == 3 && j == 3) {
    return 9;
  }
  int a = i;
  int b = j;
  if (i <= j) {
    return i * 3 + j - std::max(0, a - 1);
  }
  return j * 3 + i - std::max(0, b - 1);
}

inline bool ContainsInvalidBase(const Tuple& bases) {
  // if a pair contains an 'N'
  idx_t b1, b2, b3, b4;
  std::tie(b1, b2, b3, b4) = bases;
  if (b1 == 4 || b2 == 4 || b3 == 4 || b4 == 4) {
    return true;
  }
  return false;
}

inline idx_t CountOuterBranches(const vector<IntPair>& pairings) {
  idx_t prev_pos = -1;
  idx_t outer_branch_counter = 0;
  vector<IntPair>::const_iterator it;
  for (it = pairings.begin(); it != pairings.end(); ++it) {
    if (it->first > prev_pos) {
      outer_branch_counter++;
      prev_pos = it->second;
    }
  }
  if (outer_branch_counter == 0) {
    outer_branch_counter++;
  }
  return outer_branch_counter;
}

void InitializeFeatures(FeatureVec* vec);
void ConvertFeatures(FeatureVec* vec);

void CalculateFeatures(const Structure& structure, FeatureVec* features);

void GenericCounter(ScoringFunction_t func, const Structure& rna,
                    const Tuple& t, FeatureVec* features);

score_t GenericScorer(ScoringFunction_t func, const Structure& rna,
                      const Tuple& t, const FeatureVec& params);

#define WEIGHT_FUNC_ARGS \
  const Tuple &tp, const Structure &rna, IndexVec *indices
void HairpinWeights(WEIGHT_FUNC_ARGS);
void HelixBasePairWeights(WEIGHT_FUNC_ARGS);
void HelixChangeWeights(WEIGHT_FUNC_ARGS);
void HelixClosingWeights(WEIGHT_FUNC_ARGS);
void HelixExtendWeights(WEIGHT_FUNC_ARGS);
void HelixStackingWeights(WEIGHT_FUNC_ARGS);
void MultiBaseWeights(WEIGHT_FUNC_ARGS);
void MultiMismatchWeights(WEIGHT_FUNC_ARGS);
void MultiPairedWeights(WEIGHT_FUNC_ARGS);
void MultiUnpairedWeights(WEIGHT_FUNC_ARGS);
void OuterUnpairedWeights(WEIGHT_FUNC_ARGS);
void OuterBranchWeights(WEIGHT_FUNC_ARGS);
void SingleWeights(WEIGHT_FUNC_ARGS);

}  // namespace qrisp

#endif  // QRISP_MODEL_H_
