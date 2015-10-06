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

#ifndef QRISP_RECURRENCES_NBEST_H_
#define QRISP_RECURRENCES_NBEST_H_

#include "model.h"
#include "recurrences.h"
#include "utils.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <stack>
#include <tuple>
#include <vector>

namespace qrisp {
namespace nbest {

const int max_nbest = 2;

// A Cell corresponds to one entry in a 'classical' dynamic programming table.
// A Cell consists of
// * current state,
// * positions,
// * score, and
// * pointer to its previous cell (traceback).
// The current state and the positions are the coordinates. New Cells are
// created from old cells at either one given coordinate or two coordinates.
class Cell {
 public:
  Cell(int state, int i, int j, double s, Cell* prev)
      : state_(state),
        i_(i),
        j_(j),
        score_(s),
        trace_(prev),
        trace_snd_(nullptr) {}

  Cell(int state, int i, int j, double s, Cell* prev, Cell* prev_snd)
      : state_(state),
        i_(i),
        j_(j),
        score_(s),
        trace_(prev),
        trace_snd_(prev_snd) {}

  ~Cell() {}

  int state_;
  int i_;
  int j_;
  double score_;

  Cell* trace_;
  Cell* trace_snd_;
};

typedef std::vector<Cell*> NBestCells;

class DPTable {
 public:
  DPTable(int d1, int d2, int d3) : d1_size(d1), d2_size(d2), d3_size(d3) {
    data_.assign(d1 * d2 * d3, NBestCells());
  }

  ~DPTable();

  NBestCells& operator()(const int& i, const int& j, const int& k);
  const NBestCells& Get(const int& i, const int& j, const int& k) const;

  int GetSize1() const { return d1_size; }
  int GetSize2() const { return d2_size; }
  int GetSize3() const { return d3_size; }
  int size() const { return data_.size(); }

 private:
  int d1_size;
  int d2_size;
  int d3_size;
  vector<NBestCells> data_;
};

inline bool CellComp(const Cell* a, const Cell* b) {
  return std::tie(a->state_, a->i_, a->j_) < std::tie(b->state_, b->i_, b->j_);
}

inline bool CellCompWithScore(const Cell* a, const Cell* b) {
  return std::tie(a->state_, a->i_, a->j_, a->score_) <
         std::tie(b->state_, b->i_, b->j_, b->score_);
}

inline bool CellCompWithScoreReversed(const Cell* a, const Cell* b) {
  return std::tie(a->state_, a->i_, a->j_, a->score_) >
         std::tie(b->state_, b->i_, b->j_, b->score_);
}

inline bool operator==(const Cell& a, const Cell& b) {
  return (!CellCompWithScore(&a, &b)) && (!CellCompWithScore(&b, &a));
}

inline bool Equiv(const NBestCells& cells_a, const NBestCells& cells_b) {
  if (cells_a.size() != cells_b.size()) {
    return false;
  }
  for (int i = 0; i < cells_a.size(); i++) {
    if (CellCompWithScore(cells_a[i], cells_b[i]) ||
        CellCompWithScore(cells_b[i], cells_a[i])) {
      return false;
    }
  }
  return true;
}

void CreateCellsFromSingleCoordinate(DPTable* tbl, const int k, const int i,
                                     const int j, const int kp, const int ip,
                                     const int jp, const double score);

void CreateCellsFromTwoCoordinates(DPTable* tbl, const int k, const int i,
                                   const int j, const int kp, const int ip,
                                   const int jp, const int kpp, const int ipp,
                                   const int jpp, const double score);

void FillTable(const Structure& input, DECODING_MODE vmode, double loss_factor,
               const FeatureVec& params, DPTable* tbl);

double PerformTraceback(const int n, const DPTable& trace, const int rna_size,
                        vector<idx_t>* result);

double DecodePaths(const int n, const Structure& rna, const FeatureVec& model,
                   const DECODING_MODE& vmode, const double& lossFactor,
                   vector<vector<idx_t>>* result);

struct ObjectDeleter {
  template <typename T>
  void operator()(const T* ptr) const {
    delete ptr;
  }
};

}  // namespace nbest
}  // namespace qrisp
#endif  // QRISP_RECURRENCES_NBEST_H_
