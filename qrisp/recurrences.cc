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

#include "recurrences.h"

//#define _SHOW_STATES_
#ifdef _SHOW_STATES_
#define message(fmt, ...) printf(fmt, __VA_ARGS__)
#else
#define message(fmt, ...)
#endif

namespace qrisp {

void UpdateMax(DPTable* t1, DPTable* t2, const idx_t i, const idx_t j,
               const idx_t k, const double v1, const double v2) {
  if (v1 > (*t1)(i, j, k)) {
    (*t1)(i, j, k) = v1;
    CHECK(v2 >= 0);
    (*t2)(i, j, k) = v2;
  }
}

// Functor that hides the pairings structure and just returns loss terms for
// pairs and free residues.
// TODO(fdb): Factor out the different cases to separate functions and then use
// a function pointer defined outside the loop.
std::tuple<score_t, score_t> LossFunctor::GetLossTerms(const idx_t i,
                                                       const idx_t j) {
  double pair = 0.0;
  double single = 0.0;
  if (dmode == LOSS_ENABLED) {
    pair = (pairings_[i + 1] != 0) ? 0.0 : loss_;
    pair += (pairings_[j] != 0) ? 0.0 : loss_;
    single = (pairings_[i + 1] == 0) ? 0.0 : loss_;
  } else {
    if (dmode == PAIRWISE_LOSS) {
      pair = (pairings_[i + 1] == j) ? 0.0 : 2 * loss_;
      single = (pairings_[i + 1] == 0) ? 0.0 : loss_;
    } else {
      if (dmode == LOSS_DISABLED) {
        // TODO(fdb): why was this here in the first place?
        // pair = 1.0 * loss_;
        pair = 0.0;
        single = 0.0;
      } else {
        CHECK(false);
      }
    }
  }
  return std::make_tuple(pair, single);
}

#define BIND_SC(func_ptr) bind(GenericScorer, func_ptr, input, ph, params);
// This method calculates the entries of the viterbi and traceback tables.
void FillTables(const Structure& input, DECODING_MODE vmode, double loss_factor,
                const FeatureVec& params, DPTable* scores, DPTable* traceback) {
  // Initialize scoring functions.
  auto ph = std::placeholders::_1;
  auto Hairpin = BIND_SC(&HairpinWeights);
  auto HelixBasePair = BIND_SC(&HelixBasePairWeights);
  auto HelixChange = BIND_SC(&HelixChangeWeights);
  auto HelixClosing = BIND_SC(&HelixClosingWeights);
  auto HelixExtend = BIND_SC(&HelixExtendWeights);
  auto HelixStacking = BIND_SC(&HelixStackingWeights);
  auto MultiBase = BIND_SC(&MultiBaseWeights);
  auto MultiMismatch = BIND_SC(&MultiMismatchWeights);
  auto MultiPaired = BIND_SC(&MultiPairedWeights);
  auto MultiUnpaired = BIND_SC(&MultiUnpairedWeights);
  auto OuterBranch = BIND_SC(&OuterBranchWeights);
  auto OuterUnpaired = BIND_SC(&OuterUnpairedWeights);
  auto Single = BIND_SC(&SingleWeights);
  const idx_t rna_size = input.Size();
  const idx_t length = rna_size - 1;
  // Create the functor for losses based on the input structure.
  LossFunctor* loss_functor = new LossFunctor(input, vmode, loss_factor);
  auto UPDATE = bind(UpdateMax, scores, traceback, std::placeholders::_1,
                     std::placeholders::_2, std::placeholders::_3,
                     std::placeholders::_4, std::placeholders::_5);
  auto EncodeTB = bind(EncodeTraceback, std::placeholders::_1,
                       std::placeholders::_2, rna_size);
  // auto loss_function = bind();
  // if (i+1 > 0 && j > 0 && i+1 < rna_size && j < rna_size &&
  // input.quality(i+1) < 0.9 && input.quality(j) < 0.9) {
  //   pair += .5;
  // }
  // Some scoring functions do not need positional information. This tuple
  // provides neutral positional information for these cases.
  const Tuple ZERO = Tuple(0, 0, 0, 0);
  const idx_t D = 5;
  score_t sum = 0.0;
  for (idx_t m = 0; m <= length; m++) {
    for (idx_t i = 0, j = m; j <= length; i++, j++) {
      double pair = 0;
      double single = 0;
      std::tie(pair, single) = loss_functor->GetLossTerms(i, j);
      // DO_HELIX 0
      if (i + 2 <= j) {
        sum = HelixChange(mt3(i + 1, j, 1)) + HelixClosing(mt2(i + 1, j)) +
              pair + HelixBasePair(mt2(i + 1, j)) +
              (*scores)(DO_HELIX + 1, i + 1, j - 1);
        UPDATE(DO_HELIX, i, j, sum, EncodeTB(0, 0));
      }
      // OUTER
      if (j == length) {
        if (i == length) {
          UPDATE(DO_OUTER, i, j, 0.0, EncodeTB(0, 0));
        }
        if (i < length) {
          sum = single + OuterUnpaired(mt1(i + 1)) +
                (*scores)(DO_OUTER, i + 1, j);
          UPDATE(DO_OUTER, i, j, sum, EncodeTB(1, 0));
        }
        for (idx_t ip = i + 2; ip <= length; ip++) {
          sum = OuterBranch(ZERO) + (*scores)(DO_HELIX, i, ip) +
                (*scores)(DO_OUTER, ip, j);
          UPDATE(DO_OUTER, i, j, sum, EncodeTB(ip, j));
        }
      }
      // DO_MULTI
      if (i == j) {
        UPDATE(DO_MULTI + 2, i, j, 0.0, EncodeTB(0, 0));
      }
      for (idx_t n = 0; n <= 2; n++) {
        if (i < j) {
          sum = single + MultiUnpaired(mt3(i + 1, j, i + 1)) +
                (*scores)(DO_MULTI + n, i + 1, j);
          UPDATE(DO_MULTI + n, i, j, sum, EncodeTB(1, 0));
        }
        if (i > 0 && j < length) {
          for (idx_t jp = i + 2; jp <= j; jp++) {
            sum = MultiPaired(mt2(i, j)) +
                  MultiMismatch(mt(jp, i + 1, jp + 1, i)) +
                  (*scores)(DO_HELIX, i, jp) +
                  (*scores)(DO_MULTI + min(idx_t(2), n + 1), jp, j);
            UPDATE(DO_MULTI + n, i, j, sum, EncodeTB(jp, 0));
          }
        }
      }
      // DO_LOOP
      if (i + 3 <= j && i > 0 && j < length) {
        double hairpin_loss = loss_functor->PairingLoss(i, j);
        UPDATE(DO_LOOP, i, j, hairpin_loss + Hairpin(mt(i, j + 1, i + 1, j)),
               EncodeTB(0, 0));
        for (idx_t li = 0; li <= MAX_LOOP_SIZE; li++) {
          for (idx_t lj = 0; li + lj <= MAX_LOOP_SIZE; lj++) {
            if (li + lj > 0) {
              idx_t ip = i + li;
              idx_t jp = j - lj;
              if (!(ip + 2 <= jp)) {
                continue;
              }
              double single_loop_loss = loss_functor->PairingLoss(i, ip) +
                                        loss_functor->PairingLoss(jp, j);
              sum = single_loop_loss + Single(mt(i, j, ip, jp)) +
                    (*scores)(DO_HELIX, ip, jp);
              UPDATE(DO_LOOP, i, j, sum, EncodeTB(ip, jp));
            }
          }
        }
        if (i + 2 <= j) {
          sum = MultiBase(ZERO) + MultiPaired(ZERO) +
                MultiMismatch(mt(i, j + 1, i + 1, j)) +
                (*scores)(DO_MULTI, i, j);
          UPDATE(DO_LOOP, i, j, sum, EncodeTB(1, 1));
        }
      }

      // DO_HELIX
      if (i > 0 && j < length) {
        for (idx_t n = 1; n <= D; n++) {
          if (i + 2 <= j) {
            score_t score0 = HelixStacking(mt(i, j + 1, i + 1, j)) + pair +
                             HelixBasePair(mt2(i + 1, j));
            if (n < D) {
              sum = score0 + HelixChange(mt3(i + 1, j, n + 1)) +
                    (*scores)(DO_HELIX + n + 1, i + 1, j - 1);
              UPDATE(DO_HELIX + n, i, j, sum, EncodeTB(0, 0));
            } else {
              sum = score0 + HelixExtend(mt2(i + 1, j)) +
                    (*scores)(DO_HELIX + D, i + 1, j - 1);
              UPDATE(DO_HELIX + n, i, j, sum, EncodeTB(1, 0));
            }
          }
          if (n > 1) {
            sum = HelixClosing(mt2(j + 1, i)) + (*scores)(DO_LOOP, i, j);
            UPDATE(DO_HELIX + n, i, j, sum, EncodeTB(2, 0));
          }
        }
      }
    }
  }
  delete loss_functor;
}

void PrintTable(const DPTable& t) {
  for (idx_t s = 0; s < 8; s++) {
    for (idx_t i = 0; i < t.GetSize2(); i++) {
      for (idx_t j = 0; j < t.GetSize3(); j++) {
        if (t.Get(s, i, j) < -9999.0) {
          printf("%.03f ", -9999.0);
        } else {
          printf("%.03f ", t.Get(s, i, j));
        }
      }
      printf("\n");
    }
    printf("\n\n");
  }
}

void PerformTraceback(const DPTable& trace, const idx_t& rna_size,
                      vector<idx_t>* result) {
  // PrintTable(trace);
  result->clear();
  result->resize(rna_size, 0);
  (*result)[0] = IDX_NOT_SET;
  queue<std::tuple<idx_t, idx_t, idx_t>> path;
  // Auxiliary function that saves some code.
  auto push_triple = [&](const idx_t& x, const idx_t& y, const idx_t& z) {
    path.push(std::make_tuple(x, y, z));
  };
  push_triple(DO_OUTER, 0, rna_size - 1);
  // Enter the main loop to decode the best structure.
  while (!path.empty()) {
    idx_t v, i, j;
    std::tie(v, i, j) = path.front();
    message("Front triple: (%d, %d, %d)\n", v, i, j);
    path.pop();
    idx_t p, k;
    std::tie(p, k) = DecodeTraceback(trace.Get(v, i, j), rna_size);
    message("Traceback tuple: (%d, %d)\n", p, k);
    CHECK(i >= 0);
    CHECK(j >= 0);
    CHECK(p >= 0);
    CHECK(k >= 0);
    switch (v) {
      case DO_OUTER:
        message("DO_OUTER @ (%3d, %3d),(%1d, %1d)\n", i, j, p, k);
        if (p == 1 && k == 0) {
          push_triple(DO_OUTER, i + 1, j);
        } else if (p > 0) {
          push_triple(DO_HELIX, i, p);
          push_triple(DO_OUTER, p, j);
        }
        break;

      case DO_HELIX:
        message("DO_HELIX @ (%3d, %3d)\n", i, j);
        push_triple(DO_HELIX + 1, i + 1, j - 1);
        (*result)[i + 1] = j;
        (*result)[j] = i + 1;
        break;

      case DO_HELIX1:
      case DO_HELIX2:
      case DO_HELIX3:
      case DO_HELIX4:
      case DO_HELIX5:
        message("DO_HELIX-%d @ (%3d, %3d)\n", (v - DO_HELIX1 + 1), i, j);
        if (p == 0) {
          push_triple(v + 1, i + 1, j - 1);
          (*result)[i + 1] = j;
          (*result)[j] = i + 1;
        } else if (p == 1) {
          push_triple(v, i + 1, j - 1);
          (*result)[i + 1] = j;
          (*result)[j] = i + 1;
        } else {
          push_triple(DO_LOOP, i, j);
        }
        break;

      case DO_LOOP:
        message("DO_LOOP @ (%3d, %3d)\t", i, j);
        if (p == 0 && k == 0) {
          message("(hairpin) i,j=%d,%d \n", i, j);
          // hairpin loop
        } else if (p == 1 && k == 1) {
          message("(multi) i,j=%d,%d \n", i, j);
          push_triple(DO_MULTI, i, j);
        } else {
          message("(single) ip,jp=%d,%d \n", p, k);
          push_triple(DO_HELIX, p, k);
        }
        break;

      case DO_MULTI:
      case DO_MULTI1:
      case DO_MULTI2:
        message("DO_MULTI-%d @ (%3d, %3d)\n", (v - DO_MULTI + 1), i, j);
        if (p == 1 && k == 0) {
          push_triple(v, i + 1, j);
        } else if (p >= 2) {
          push_triple(DO_HELIX, i, p);
          push_triple(DO_MULTI + min(idx_t(2), v - DO_MULTI + 1), p, j);
        }
        break;

      default:
        CHECK(false);
    };
  }
}

double DecodeBestPath(const Structure& rna, const FeatureVec& model,
                      const DECODING_MODE& vmode, const double& loss_factor,
                      vector<idx_t>* result) {
  const idx_t rna_size = rna.Size();
  FeatureVec model_transformed = model;
  InitializeFeatures(&model_transformed);
  // Create and fill dynamic programming tables.
  DPTable* scores = new DPTable(NUM_STATES, rna_size, rna_size, LOG_ZERO);
  DPTable* traceback = new DPTable(NUM_STATES, rna_size, rna_size, LOG_ZERO);
  FillTables(rna, vmode, loss_factor, model_transformed, scores, traceback);
  // Perform backtracking on the tables to infer the highest scoring path.
  PerformTraceback(*traceback, rna_size, result);
  double best_path_score = scores->Get(DO_OUTER, 0, rna_size - 1);
  delete scores;
  delete traceback;
  return best_path_score;
}

}  // namespace qrisp
