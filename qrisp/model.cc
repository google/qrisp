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

#include <functional>
#include <set>

#include "model.h"
#include "plif.h"
#include "qrisp/proto/parameters.pb.h"
#include "third_party/googleflags/include/gflags/gflags.h"
#include "utils.h"

//#define _SHOW_STATES_
#ifdef _SHOW_STATES_
#define message(fmt, ...) printf(fmt, __VA_ARGS__)
#else
#define message(fmt, ...)
#endif

static constexpr bool enable_quality_features = false;

namespace qrisp {

void UpdateFeatureWeights(const IndexVec& updates, FeatureVec* score_vec) {
  for (auto it = updates.begin(); it != updates.end(); ++it) {
    auto insert_result = score_vec->insert(*it);
    if (!insert_result.second) {
      (*score_vec)[it->first] += it->second;
    }
  }
}

void GenericCounter(ScoringFunction_t func, const Structure& rna,
                    const Tuple& t, FeatureVec* features) {
  IndexVec indices;
  (*func)(t, rna, &indices);
  UpdateFeatureWeights(indices, features);
}

score_t GetSummedParameters(const IndexVec& counts, const FeatureVec& scores) {
  score_t result = 0.0;
  for (auto it = counts.begin(); it != counts.end(); ++it) {
    if (scores.find(it->first) != scores.end()) {
      result += scores.at(it->first) * it->second;
    }
  }
  return result;
}

score_t GenericScorer(ScoringFunction_t func, const Structure& rna,
                      const Tuple& t, const FeatureVec& params) {
  IndexVec indices;
  (*func)(t, rna, &indices);
  return GetSummedParameters(indices, params);
}

inline void push(IndexVec* indices, idx_t i, score_t s) {
  indices->push_back(make_pair(i, s));
}

inline void push_one(IndexVec* indices, idx_t i) {
  indices->push_back(make_pair(i, 1.0));
}

void HairpinWeights(const Tuple& idx, const Structure& rna, IndexVec* indices) {
  const idx_t hp_len = std::get<3>(idx) - std::get<0>(idx);
  const Tuple& bases = rna.GetBasesAt(idx);
  push_one(indices,
           Features::TERMINAL_MISMATCH_AAAA + CalculateFullOffset(bases));
  if (hp_len <= 30) {
    push_one(indices, Features::HAIRPIN_LENGTH_0 + hp_len);
  } else {
    push_one(indices, Features::HAIRPIN_LENGTH_30);
    push(indices, Features::HAIRPIN_EXTEND, hp_len - 30.0);
  }
  if (enable_quality_features) {
    for (int k = std::get<2>(idx); k <= std::get<3>(idx); k++) {
      idx_t feature_offset =
          Features::PLIF1D_OFFSET + rna[k] * Features::PLIF_INCREMENT;
      EvaluatePlifAt(limits, rna.quality(k), feature_offset, indices);
    }
  }
}

// TODO(fabio): Add a check after GetBasesAt that makes sure only bases that
// were set are used.
void HelixBasePairWeights(const Tuple& idx, const Structure& rna,
                          IndexVec* indices) {
  const Tuple& bases = rna.GetBasesAt(idx);
  const idx_t offset =
      GetBasePairOffset(std::get<0>(bases), std::get<1>(bases));
  push_one(indices, Features::BASE_PAIR_AA + offset);
  if (enable_quality_features) {
    set<idx_t> possible_offsets = {2000, 2040, 2080, 2120, 2160,
                                   2200, 2240, 2280, 2320, 2360};
    idx_t feature_offset =
        Features::PLIF2D_OFFSET + offset * Features::PLIF_INCREMENT;
    assert(possible_offsets.find(feature_offset) != possible_offsets.end());
    Evaluate2DPlifAt(limits, limits, rna.quality(std::get<0>(idx)),
                     rna.quality(std::get<1>(idx)), feature_offset, indices);
    // idx_t feature_offset = PLIF2D_OFFSET + rna[get<0>(idx)] * PLIF_INCREMENT;
    // EvaluatePlifAt(limits, rna.quality(get<0>(idx)), feature_offset,
    // indices);
    ////
    // feature_offset = PLIF2D_OFFSET + rna[get<1>(idx)] * PLIF_INCREMENT;
    // EvaluatePlifAt(limits, rna.quality(get<1>(idx)), feature_offset,
    // indices);
  }
}

void HelixChangeWeights(const Tuple& idx, const Structure& rna,
                        IndexVec* indices) {
  const idx_t len = std::get<2>(idx);
  assert(len > 0);
  if (len <= 5) {
    push_one(indices, Features::STACKING_1_SCORE + len - 1);
  } else {
    push_one(indices, Features::EXTENSION_SCORE);
  }
}

void HelixClosingWeights(const Tuple& idx, const Structure& rna,
                         IndexVec* indices) {
  const Tuple& bases = rna.GetBasesAt(idx);
  CHECK(std::get<0>(bases) != 100);
  CHECK(std::get<1>(bases) != 100);
  const idx_t offset = std::get<0>(bases) * 4 + std::get<1>(bases);
  push_one(indices, Features::HELIX_CLOSING_AA + offset);
}

void HelixExtendWeights(const Tuple& idx, const Structure& rna,
                        IndexVec* indices) {
  push_one(indices, Features::EXTENSION_SCORE);
}

void HelixStackingWeights(const Tuple& idx, const Structure& rna,
                          IndexVec* indices) {
  const Tuple& bases = rna.GetBasesAt(idx);
  push_one(indices, Features::STACKING_PAIR_AAAA + CalculateOffset(bases));
}

void MultiBaseWeights(const Tuple& idx, const Structure& rna,
                      IndexVec* indices) {
  push_one(indices, Features::MULTI_BASE);
}

void MultiMismatchWeights(const Tuple& idx, const Structure& rna,
                          IndexVec* indices) {
  const Tuple& bases = rna.GetBasesAt(idx);
  push_one(indices, Features::SINGLE_BASE_PAIR_STACKING_LEFT_AAA +
                        CalculateLeftOffset(bases));
  push_one(indices, Features::SINGLE_BASE_PAIR_STACKING_RIGHT_AAA +
                        CalculateRightOffset(bases));
}

void MultiPairedWeights(const Tuple& idx, const Structure& rna,
                        IndexVec* indices) {
  push_one(indices, Features::MULTI_PAIRED);
}

void MultiUnpairedWeights(const Tuple& idx, const Structure& rna,
                          IndexVec* indices) {
  push_one(indices, Features::MULTI_UNPAIRED);
  if (enable_quality_features) {
    const Tuple bases = rna.GetBasesAt(idx);
    const idx_t feature_offset =
        Features::PLIF1D_OFFSET + std::get<0>(bases) * Features::PLIF_INCREMENT;
    EvaluatePlifAt(limits, rna.quality(std::get<0>(idx)), feature_offset,
                   indices);
  }
}

void OuterUnpairedWeights(const Tuple& idx, const Structure& rna,
                          IndexVec* indices) {
  push_one(indices, Features::OUTER_UNPAIRED);
  if (enable_quality_features) {
    const Tuple bases = rna.GetBasesAt(idx);
    const idx_t feature_offset =
        Features::PLIF1D_OFFSET + std::get<0>(bases) * Features::PLIF_INCREMENT;
    EvaluatePlifAt(limits, rna.quality(std::get<0>(idx)), feature_offset,
                   indices);
  }
}

void OuterBranchWeights(const Tuple& idx, const Structure& rna,
                        IndexVec* indices) {
  push_one(indices, Features::OUTER_BRANCH);
}

void SingleWeights(const Tuple& idx, const Structure& rna, IndexVec* indices) {
  idx_t i, j, ip, jp;
  std::tie(i, j, ip, jp) = idx;
  // assert(i > 0);
  // assert(ip > 0);
  // assert(j > 0);
  // assert(jp > 0);
  // const idx_t rna_size = rna.Size();
  // assert(i < rna_size);
  // assert(ip < rna_size);
  // assert(j < rna_size);
  // assert(jp < rna_size);
  // Check soundness of pairs
  // assert(i <= ip);
  // assert( jp <= j);
  const idx_t downstream_size = ip - i;
  const idx_t upstream_size = j - jp;
  const idx_t loop_size = downstream_size + upstream_size;
  assert(loop_size <= 30 && loop_size > 0);
  const Tuple& bases1 = rna.GetBasesAt(Tuple(i, j + 1, i + 1, j));
  push_one(indices,
           Features::TERMINAL_MISMATCH_AAAA + CalculateFullOffset(bases1));
  const Tuple& bases2 = rna.GetBasesAt(Tuple(ip + 1, jp, ip, jp + 1));
  push_one(indices,
           Features::TERMINAL_MISMATCH_AAAA + CalculateFullOffset(bases2));
  //
  if (enable_quality_features) {
    for (int k = i + 1; k <= ip; k++) {
      idx_t feature_offset =
          Features::PLIF1D_OFFSET + rna[k] * Features::PLIF_INCREMENT;
      EvaluatePlifAt(limits, rna.quality(k), feature_offset, indices);
    }
    for (int k = jp + 1; k <= j; k++) {
      // Tuple bases = rna.GetBasesAt(mt1(k));
      idx_t feature_offset =
          Features::PLIF1D_OFFSET + rna[k] * Features::PLIF_INCREMENT;
      EvaluatePlifAt(limits, rna.quality(k), feature_offset, indices);
    }
  }
  // Check whether we have a bulge loop
  if (j - jp == 0 || ip - i == 0) {
    push_one(indices, Features::BULGE_LENGTH_0 + loop_size);
  } else {
    push_one(indices, Features::INTERNAL_LENGTH_1 + loop_size - 1);
    int asym_small = downstream_size;
    int asym_big = upstream_size;
    if (downstream_size > upstream_size) {
      asym_small = upstream_size;
      asym_big = downstream_size;
    }
    int asym_offset = 0;
    for (idx_t r = 1; r < asym_small; r++) {
      asym_offset += (30 - r) - r;
    }
    push_one(indices,
             Features::INTERNAL_FULL_1_1 + (asym_offset + asym_big) - 1);
    push_one(indices,
             Features::INTERNAL_ASYMMETRY_0 + std::abs(asym_big - asym_small));
  }
}

void InitializeFeatures(FeatureVec* vec) {
  for (int i = 1; i <= MAX_LOOP_SIZE; i++) {
    (*vec)[Features::HAIRPIN_LENGTH_0 + i] +=
        (*vec)[Features::HAIRPIN_LENGTH_0 + i - 1];
    (*vec)[Features::BULGE_LENGTH_0 + i] +=
        (*vec)[Features::BULGE_LENGTH_0 + i - 1];
    (*vec)[Features::INTERNAL_ASYMMETRY_0 + i] +=
        (*vec)[Features::INTERNAL_ASYMMETRY_0 + i - 1];
    if (i < MAX_LOOP_SIZE) {
      (*vec)[Features::INTERNAL_LENGTH_1 + i] +=
          (*vec)[Features::INTERNAL_LENGTH_1 + i - 1];
    }
  }
}

void ConvertFeatures(FeatureVec* vec) {
  for (int i = 1; i < MAX_LOOP_SIZE; i++) {
    for (int j = 0; j < i; j++) {
      (*vec)[Features::HAIRPIN_LENGTH_0 + j] +=
          (*vec)[Features::HAIRPIN_LENGTH_0 + i];
      (*vec)[Features::BULGE_LENGTH_0 + j] +=
          (*vec)[Features::BULGE_LENGTH_0 + i];
      (*vec)[Features::INTERNAL_ASYMMETRY_0 + j] +=
          (*vec)[Features::INTERNAL_ASYMMETRY_0 + i];
      (*vec)[Features::INTERNAL_LENGTH_1 + j] +=
          (*vec)[Features::INTERNAL_LENGTH_1 + i];
    }
  }
  for (int j = 0; j < MAX_LOOP_SIZE; j++) {
    (*vec)[Features::HAIRPIN_LENGTH_0 + j] +=
        (*vec)[Features::HAIRPIN_LENGTH_0 + MAX_LOOP_SIZE];
    (*vec)[Features::BULGE_LENGTH_0 + j] +=
        (*vec)[Features::BULGE_LENGTH_0 + MAX_LOOP_SIZE];
    (*vec)[Features::INTERNAL_ASYMMETRY_0 + j] +=
        (*vec)[Features::INTERNAL_ASYMMETRY_0 + MAX_LOOP_SIZE];
  }
}

#define BIND(func_ptr) bind(GenericCounter, func_ptr, structure, ph, features)
void CalculateFeatures(const Structure& structure, FeatureVec* features) {
  features->clear();
  // Initialize scoring functions
  auto ph = std::placeholders::_1;
  auto HairpinCounter = BIND(&HairpinWeights);
  auto HelixBasePairCounter = BIND(&HelixBasePairWeights);
  auto HelixChangeCounter = BIND(&HelixChangeWeights);
  auto HelixClosingCounter = BIND(&HelixClosingWeights);
  auto HelixExtendCounter = BIND(&HelixExtendWeights);
  auto HelixStackingCounter = BIND(&HelixStackingWeights);
  auto MultiBaseCounter = BIND(&MultiBaseWeights);
  auto MultiMismatchCounter = BIND(&MultiMismatchWeights);
  auto MultiPairedCounter = BIND(&MultiPairedWeights);
  auto MultiUnpairedCounter = BIND(&MultiUnpairedWeights);
  auto OuterBranchCounter = BIND(&OuterBranchWeights);
  auto OuterUnpairedCounter = BIND(&OuterUnpairedWeights);
  auto SingleCounter = BIND(&SingleWeights);

  // Now calculate the substructures.
  vector<Substructure> substructures;
  vector<bool> acc_pos;
  structure.CalculateSubstructure(&substructures, &acc_pos);
  // const int length = structure.etSize() - 1;
  // Some counting function do not need positional information. This tuple
  // provides neutral positional information for these cases.
  const Tuple zero_tuple = Tuple(1, 1, 1, 1);
  // If there is no base-pair set all pos to unpaired and return.
  message("CalculateFeatures: Substructure size: " << substructures.size()
                                                   << endl);
  if (substructures.size() == 0) {
    message("No substructures present %c", '\n');
    for (int i = 1; i < structure.Size(); i++) {
      OuterUnpairedCounter(mt1(i));
    }
    // CountOuterUnpairedInterval(1, length);
    return;
  }
  // Before we loop over the base-pairs we calculate the accessible positions.
  // Important: Skip the 0th position here. It is no used in the structure.
  int idx = 0;
  for (auto it = acc_pos.cbegin() + 1; it != acc_pos.cend(); ++it) {
    if (!*it) {
      message("outer unpaired at %d", idx);
      OuterUnpairedCounter(mt1(idx + 1));
    }
    idx++;
  }
  idx_t stacking_depth = 1;
  idx_t last_covered_position = 0;
  // Now loop over the pairings.
  for (auto jt = substructures.cbegin(); jt != substructures.cend(); ++jt) {
    idx_t num_inner_pairs = jt->GetNumberOfInnerPairs();
    idx_t i = jt->GetOuterPair().first;
    idx_t j = jt->GetOuterPair().second;
    message("i, j / num_inner_pairs: %d, %d / %d\n", i, j, num_inner_pairs);
    assert(i > 0);
    // assert(j <= length);
    // assert (i + 2 <= j);
    message("i / pos are %d,%d\n", i, last_covered_position);
    if (i > last_covered_position) {
      message("outer branch at %d,%d\n", i, j);
      OuterBranchCounter(mt2(i, j));
      last_covered_position = j;
    }
    // if (i + 2 <= j) {
    //  continue;
    //}
    // HAIRPIN
    if (num_inner_pairs == 0) {
      message("hairpin at %d,%d / depth: %d\n", i, j, stacking_depth);
      HairpinCounter(mt(i, j, i + 1, j - 1));
      HelixBasePairCounter(mt2(i, j));
      HelixChangeCounter(mt3(i, j, stacking_depth));
      HelixClosingCounter(mt2(j, i));
      stacking_depth = 1;
    }
    // STEM
    if (num_inner_pairs == 1) {
      idx_t ip = jt->GetInnerPairAt(0).first;
      idx_t jp = jt->GetInnerPairAt(0).second;
      if ((ip - 1 - i) + (j - 1 - jp) <= MAX_LOOP_SIZE) {
        //
        if (ip - i == 1 && j - jp == 1) {
          message("stem at %d,%d\t", i, j);
          HelixStackingCounter(mt(i, j, ip, jp));
          HelixBasePairCounter(mt2(i, j));
          HelixChangeCounter(mt3(i, j, stacking_depth));
          message("depth is %d\n", stacking_depth);
          if (stacking_depth == 1) {
            HelixClosingCounter(mt2(i, j));
          }
          message("Stacking at %d,%d,%d,%d\n", i - 1, j + 1, i, j);
          stacking_depth++;
        } else {
          message("interior/bulge at %d,%d %d,%d\n", i, j, ip, jp);
          SingleCounter(mt(i, j - 1, ip - 1, jp));
          HelixBasePairCounter(mt2(i, j));
          HelixChangeCounter(mt3(i, j, stacking_depth));
          HelixClosingCounter(mt2(j, i));
          stacking_depth = 1;
        }
      }
    }
    // MULTI LOOP
    if (num_inner_pairs >= 2) {
      message("multi at %d,%d / depth: %d\n", i, j, stacking_depth);
      if (stacking_depth == 1) {
        HelixClosingCounter(mt2(i, j));
      }
      HelixBasePairCounter(mt2(i, j));
      HelixChangeCounter(mt3(i, j, stacking_depth));
      HelixClosingCounter(mt2(j, i));
      stacking_depth = 1;
      MultiBaseCounter(zero_tuple);
      MultiPairedCounter(zero_tuple);
      MultiMismatchCounter(mt(i, j, i + 1, j - 1));
      idx_t closed = 0;
      idx_t outer_start = i;
      idx_t outer_stop = j;
      idx_t unpaired2 = 0;
      for (int k = 0; k < num_inner_pairs; k++) {
        idx_t ip = jt->GetInnerPairAt(k).first;
        idx_t jp = jt->GetInnerPairAt(k).second;
        outer_stop = ip;
        MultiPairedCounter(zero_tuple);
        MultiMismatchCounter(mt(jp, ip, jp + 1, ip - 1));
        closed += (jp - ip) + 1;
        for (int ii = outer_start + 1; ii < outer_stop; ii++) {
          unpaired2++;
          MultiUnpairedCounter(mt1(ii));
        }
        outer_start = jp;
      }
      outer_stop = j;
      for (int ii = outer_start + 1; ii < outer_stop; ii++) {
        unpaired2++;
        MultiUnpairedCounter(mt1(ii));
      }
      const idx_t unpaired = (j - i) - closed - 1;
      CHECK(unpaired == unpaired2);
    }
  }
  ConvertFeatures(features);
}

}  // namespace qrisp
