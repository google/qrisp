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

#include "recurrences-nbest.h"

namespace qrisp {
namespace nbest {

const NBestCells& DPTable::Get(const int& i, const int& j, const int& k) const {
  CHECK(i < d1_size && j < d2_size && k < d3_size);
  return data_[i + d1_size * (j + d2_size * k)];
}

NBestCells& DPTable::operator()(const int& i, const int& j, const int& k) {
  CHECK(i < d1_size && j < d2_size && k < d3_size);
  return data_[i + d1_size * (j + d2_size * k)];
}

DPTable::~DPTable() {
  for (auto& cells : data_) {
    for_each(cells.begin(), cells.end(), nbest::ObjectDeleter());
  }
}

void CreateCellsFromSingleCoordinate(DPTable* tbl, const int k, const int i,
                                     const int j, const int kp, const int ip,
                                     const int jp, const double score) {
  NBestCells& current_cells = (*tbl)(k, i, j);
  if (kp == -1) {
    current_cells.push_back(new Cell(k, i, j, score, nullptr));
  } else {
    // Find source cells.
    const NBestCells& previous_cells = tbl->Get(kp, ip, jp);
    // Create new cells.
    for (auto it = previous_cells.begin(); it != previous_cells.end(); ++it) {
      current_cells.push_back(new Cell(k, i, j, score + (*it)->score_, *it));
    }
  }
  if (current_cells.size() > max_nbest) {
    sort(current_cells.begin(), current_cells.end(), CellCompWithScoreReversed);
    for (int i = max_nbest; i < current_cells.size(); i++) {
      delete current_cells[i];
    }
    current_cells.resize(max_nbest);
  }
}

void CreateCellsFromTwoCoordinates(DPTable* tbl, const int k, const int i,
                                   const int j, const int kp, const int ip,
                                   const int jp, const int kpp, const int ipp,
                                   const int jpp, const double score) {
  const NBestCells& cells_a = tbl->Get(kp, ip, jp);
  const NBestCells& cells_b = tbl->Get(kpp, ipp, jpp);
  NBestCells& current_cells = (*tbl)(k, i, j);
  for (auto it = cells_a.begin(); it != cells_a.end(); ++it) {
    for (auto jt = cells_b.begin(); jt != cells_b.end(); ++jt) {
      current_cells.push_back(
          new Cell(k, i, j, (*it)->score_ + (*jt)->score_ + score, *it, *jt));
    }
  }
  if (current_cells.size() > max_nbest) {
    sort(current_cells.begin(), current_cells.end(), CellCompWithScoreReversed);
    for (int i = max_nbest; i < current_cells.size(); i++) {
      delete current_cells[i];
    }
    current_cells.resize(max_nbest);
  }
}

namespace sp = std::placeholders;

#define BIND_SC(func_ptr) bind(GenericScorer, func_ptr, input, ph, params);
void FillTable(const Structure& input, DECODING_MODE vmode, double loss_factor,
               const FeatureVec& params, DPTable* tbl) {
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
  const size_t length = rna_size - 1;
  // Create the functor for losses based on the input structure.
  LossFunctor* loss_functor = new LossFunctor(input, vmode, loss_factor);
  // Insertion functions.
  auto INSERT = std::bind(CreateCellsFromSingleCoordinate, tbl, sp::_1, sp::_2,
                          sp::_3, sp::_4, sp::_5, sp::_6, sp::_7);

  auto INSERT2 =
      std::bind(CreateCellsFromTwoCoordinates, tbl, sp::_1, sp::_2, sp::_3,
                sp::_4, sp::_5, sp::_6, sp::_7, sp::_8, sp::_9, sp::_10);
  // Some scoring functions do not need positional information. This tuple
  // provides neutral positional information for these cases.
  const Tuple ZERO = Tuple(0, 0, 0, 0);
  const int D = 5;
  score_t sum = 0.0;
  for (size_t m = 0; m <= length; m++) {
    // if (tbl->size() > 10000) {
    //  exit(0);
    //}
    // LOG(INFO) << "Outer index: " << m;
    for (size_t i = 0, j = m; j <= length; i++, j++) {
      double pair = 0;
      double single = 0;
      std::tie(pair, single) = loss_functor->GetLossTerms(i, j);
      // DO_HELIX 0
      if (i + 2 <= j) {
        sum = HelixChange(mt3(i + 1, j, 1)) + HelixClosing(mt2(i + 1, j)) +
              pair + HelixBasePair(mt2(i + 1, j));
        INSERT(DO_HELIX, i, j, DO_HELIX + 1, i + 1, j - 1, sum);
      }
      // OUTER
      if (j == length) {
        if (i == length) {
          INSERT(DO_OUTER, i, j, -1, 0, 0, 0.0);
        }
        if (i < length) {
          sum = single + OuterUnpaired(mt1(i + 1));
          INSERT(DO_OUTER, i, j, DO_OUTER, i + 1, j, sum);
        }
        for (size_t ip = i + 2; ip <= length; ip++) {
          sum = OuterBranch(ZERO);
          INSERT2(DO_OUTER, i, j, DO_HELIX, i, ip, DO_OUTER, ip, j, sum);
        }
      }
      // DO_MULTI
      if (i == j) {
        INSERT(DO_MULTI + 2, i, j, -1, 0, 0, 0.0);
      }
      for (int n = 0; n <= 2; n++) {
        if (i < j) {
          sum = single + MultiUnpaired(mt3(i + 1, j, i + 1));
          INSERT(DO_MULTI + n, i, j, DO_MULTI + n, i + 1, j, sum);
        }
        if (i > 0 && j < length) {
          for (size_t jp = i + 2; jp <= j; jp++) {
            sum = MultiPaired(mt2(i, j)) +
                  MultiMismatch(mt(jp, i + 1, jp + 1, i));
            INSERT2(DO_MULTI + n, i, j, DO_HELIX, i, jp,
                    DO_MULTI + min(2, n + 1), jp, j, sum);
          }
        }
      }
      // DO_LOOP
      if (i + 3 <= j && i > 0 && j < length) {
        double hairpin_loss = loss_functor->PairingLoss(i, j);
        INSERT(DO_LOOP, i, j, -1, 0, 0,
               hairpin_loss + Hairpin(mt(i, j + 1, i + 1, j)));
        for (int li = 0; li <= MAX_LOOP_SIZE; li++) {
          for (int lj = 0; li + lj <= MAX_LOOP_SIZE; lj++) {
            if (li + lj > 0) {
              int ip = i + li;
              int jp = j - lj;
              if (!(ip + 2 <= jp)) {
                continue;
              }
              double single_loop_loss = loss_functor->PairingLoss(i, ip) +
                                        loss_functor->PairingLoss(jp, j);
              sum = single_loop_loss + Single(mt(i, j, ip, jp));
              INSERT(DO_LOOP, i, j, DO_HELIX, ip, jp, sum);
            }
          }
        }
        if (i + 2 <= j) {
          sum = MultiBase(ZERO) + MultiPaired(ZERO) +
                MultiMismatch(mt(i, j + 1, i + 1, j));
          INSERT(DO_LOOP, i, j, DO_MULTI, i, j, sum);
        }
      }

      // DO_HELIX
      if (i > 0 && j < length) {
        for (int n = 1; n <= D; n++) {
          if (i + 2 <= j) {
            score_t score0 = HelixStacking(mt(i, j + 1, i + 1, j)) + pair +
                             HelixBasePair(mt2(i + 1, j));
            if (n < D) {
              sum = score0 + HelixChange(mt3(i + 1, j, n + 1));
              INSERT(DO_HELIX + n, i, j, DO_HELIX + n + 1, i + 1, j - 1, sum);
            } else {
              sum = score0 + HelixExtend(mt2(i + 1, j));
              INSERT(DO_HELIX + n, i, j, DO_HELIX + D, i + 1, j - 1, sum);
            }
          }
          if (n > 1) {
            sum = HelixClosing(mt2(j + 1, i));
            INSERT(DO_HELIX + n, i, j, DO_LOOP, i, j, sum);
          }
        }
      }
    }
  }
  delete loss_functor;
}

score_t PerformTraceback(const int n, const DPTable& tbl, const int rna_size,
                         vector<idx_t>* result) {
  result->clear();
  result->resize(rna_size, 0);
  (*result)[0] = IDX_NOT_SET;
  NBestCells cells(tbl.Get(DO_OUTER, 0, rna_size - 1));
  std::sort(cells.begin(), cells.end(), CellCompWithScoreReversed);
  if (n >= cells.size()) {
    return 0.0;
  }
  Cell* nth_final_cell = *(cells.begin() + n);
  if (nth_final_cell == nullptr) {
    return 0.0;
  }
  score_t score = nth_final_cell->score_;
  // Initialize stack by pushing n-best end cell.
  queue<Cell*> path;
  path.push(nth_final_cell);
  // Enter the main loop to decode the best structure.
  while (!path.empty()) {
    Cell* c = path.front();
    const int v = c->state_;
    const int i = c->i_;
    const int j = c->j_;
    path.pop();
    if (c->trace_ != nullptr) {
      path.push(c->trace_);
    }
    if (c->trace_snd_ != nullptr) {
      path.push(c->trace_snd_);
    }
    CHECK(i >= 0);
    CHECK(j >= 0);
    switch (v) {
      case DO_OUTER:
      case DO_LOOP:
      case DO_MULTI:
      case DO_MULTI1:
      case DO_MULTI2:
        break;

      case DO_HELIX:
      case DO_HELIX1:
      case DO_HELIX2:
      case DO_HELIX3:
      case DO_HELIX4:
      case DO_HELIX5:
        if (j - i > 3 && c->trace_->state_ != DO_LOOP) {
          (*result)[i + 1] = j;
          (*result)[j] = i + 1;
        }
        break;

      default:
        CHECK(false);
    };
  }
  return score;
}

double DecodePaths(const int n, const Structure& rna, const FeatureVec& model,
                   const DECODING_MODE& vmode, const double& loss_factor,
                   vector<vector<idx_t>>* result) {
  const idx_t rna_size = rna.Size();
  FeatureVec model_transformed = model;
  InitializeFeatures(&model_transformed);
  // Create and fill dynamic programming tables.
  nbest::DPTable tbl(NUM_STATES, rna_size, rna_size);
  FillTable(rna, vmode, loss_factor, model_transformed, &tbl);
  // Perform backtracking on the tables to infer the highest scoring path.
  for (int i = 0; i < n; i++) {
    vector<idx_t> res;
    PerformTraceback(i, tbl, rna_size, &res);
    result->push_back(res);
  }
  double best_path_score = 0.0;
  return best_path_score;
}

}  // namespace nbest
}  // namespace qrisp
