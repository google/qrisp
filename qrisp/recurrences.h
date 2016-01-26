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

#ifndef QRISP_RECURRENCES_H_
#define QRISP_RECURRENCES_H_

#include "model.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <queue>
#include <stack>

namespace qrisp {

// There are three modes for the dynamic program:
// - one without any loss terms used for prediction, and
// - one with loss terms added at the states for creating the worst margin
//   violators, and
// - one with pairwise loss enabled.
enum DECODING_MODE { LOSS_DISABLED = 0, LOSS_ENABLED = 1, PAIRWISE_LOSS = 3 };

// The following state names were adopted from the Contrafold project.
enum STATES {
  DO_OUTER = 0,
  DO_HELIX = 1,
  DO_HELIX1 = 2,
  DO_HELIX2 = 3,
  DO_HELIX3 = 4,
  DO_HELIX4 = 5,
  DO_HELIX5 = 6,
  DO_LOOP = 7,
  DO_MULTI = 8,
  DO_MULTI1 = 9,
  DO_MULTI2 = 10,
  NUM_STATES = 11
};

void DumpTables(const idx_t& d1, const idx_t& d2, const idx_t& d3, DPTable& t1,
                DPTable& t2, DPTable& t3);

inline score_t EncodeTraceback(const idx_t& i, const idx_t& j,
                               const idx_t& rna_size) {
  CHECK(i >= 0);
  CHECK(j >= 0);
  CHECK(rna_size > 0);
  CHECK(i <= rna_size);
  CHECK(j <= rna_size);
  return (i + 1) * (rna_size + 1) + (j + 1);
}

inline std::tuple<idx_t, idx_t> DecodeTraceback(const score_t& s,
                                                const idx_t& rna_size) {
  CHECK(s >= 0.0);
  CHECK(rna_size > 0);
  const idx_t val = (idx_t)s;
  return std::make_tuple(val / (rna_size + 1) - 1, val % (rna_size + 1) - 1);
}

void PerformTraceback(const DPTable& trace, const idx_t& rna_size,
                      vector<idx_t>* result);

double DecodeBestPath(const Structure& rna, const FeatureVec& model,
                      const DECODING_MODE& vmode, const double& lossFactor,
                      vector<idx_t>* result);

// Functor that hides the pairings structure and just returns loss terms for
// pairs and free residues.
class LossFunctor {
 public:
  LossFunctor(const Structure& structure, DECODING_MODE dm, score_t newloss)
      : dmode(dm), loss_(newloss) {
    structure.GetPairings(&pairings_);
    for (size_t i = 0; i < structure.Size(); i++) {
      for (size_t j = i + 1; j < structure.Size(); j++) {
        double loss_term = 0.0;
        for (size_t z = i; z < j; z++) {
          if (pairings_[z + 1] != 0) {
            loss_term += loss_;
          }
        }
        pairings_loss_.insert(make_pair(make_pair(i, j), loss_term));
      }
    }
  }

  // This method calculates the loss terms for pairs that occur withing the
  // half-open interval [i,j).
  std::tuple<score_t, score_t> GetLossTerms(const idx_t i, const idx_t j);

  // This method calculates the loss terms for pairs that occur withing the
  // half-open interval [i,j).
  inline double PairingLoss(idx_t i, idx_t j) {
    if (dmode == LOSS_DISABLED) {
      return 0.0;
    }
    return pairings_loss_[make_pair(i, j)];
  }

 private:
  DECODING_MODE dmode;
  vector<idx_t> pairings_;
  map<pair<idx_t, idx_t>, idx_t> pairings_loss_;
  score_t loss_;
};

void FillTables(const Structure& input, DECODING_MODE vmode, double loss_factor,
                const FeatureVec& params, DPTable* scores, DPTable* traceback);

void UpdateMax(DPTable* t1, DPTable* t2, const idx_t i, const idx_t j,
               const idx_t k, const double v1, const double v2);

}  // namespace qrisp

#endif  // QRISP_RECURRENCES_H_
