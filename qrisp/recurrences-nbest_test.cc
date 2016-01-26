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

#include "experimental/users/fdb/qrsp/recurrences-nbest.h"
#include "experimental/users/fdb/qrsp/proto/parameters.pb.h"
#include "experimental/users/fdb/qrsp/recurrences.h"
#include "experimental/users/fdb/qrsp/utils.h"
#include "testing/base/public/gmock.h"
#include "testing/base/public/googletest.h"
#include "testing/base/public/gunit.h"

namespace fdb {
namespace qrsp {
namespace {

using nbest::Cell;
// See http://goto/gunitprimer for an introduction to gUnit.

class RecurrencesNbestTest : public ::testing::Test {
 protected:
  RecurrencesNbestTest() {}
  ~RecurrencesNbestTest() override {}
};

TEST_F(RecurrencesNbestTest, CellComparison) {
  auto* c1 = new Cell(1, 2, 3, 0.12, nullptr);
  auto* c2 = new Cell(1, 2, 3, 0.53, nullptr);
  auto* c3 = new Cell(10, 2, 3, 0.71, nullptr);
  auto* c4 = new Cell(1, 2, 33, 0.71, nullptr);

  EXPECT_FALSE(CellComp(c1, c1));
  EXPECT_FALSE(CellCompWithScore(c1, c1));

  EXPECT_FALSE(CellComp(c1, c2));
  EXPECT_TRUE(CellCompWithScore(c1, c2));

  EXPECT_TRUE(CellComp(c2, c3));
  EXPECT_TRUE(CellCompWithScore(c2, c3));

  EXPECT_TRUE(CellComp(c1, c3));
  EXPECT_TRUE(CellCompWithScore(c1, c3));

  EXPECT_TRUE(CellComp(c1, c4));
  EXPECT_TRUE(CellCompWithScore(c1, c4));

  EXPECT_FALSE(CellComp(c4, c1));
  EXPECT_FALSE(CellCompWithScore(c4, c1));

  EXPECT_EQ(c1, c1);

  delete c1;
  delete c2;
  delete c3;
  delete c4;
}

// Test the single cell insertion. A new cell can be reached from n previous
// cells and lead to n new cells.
TEST_F(RecurrencesNbestTest, CellInserter) {
  const int rna_size = 40;
  int k, i, j;
  std::tie(k, i, j) = std::make_tuple(1, 2, 3);
  int kp, ip, jp;
  std::tie(kp, ip, jp) = std::make_tuple(9, 22, 33);

  nbest::DPTable tbl(NUM_STATES, rna_size, rna_size);
  tbl(k, i, j).push_back(new Cell(k, i, j, 0.12, nullptr));
  tbl(k, i, j).push_back(new Cell(k, i, j, 0.53, nullptr));
  tbl(k, i, j).push_back(new Cell(k, i, j, 0.71, nullptr));
  tbl(k, i, j).push_back(new Cell(k, i, j, 0.88, nullptr));
  nbest::CreateCellsFromSingleCoordinate(&tbl, kp, ip, jp, k, i, j, 0.77);

  nbest::DPTable expected_tbl(NUM_STATES, rna_size, rna_size);
  auto* c1 = new Cell(k, i, j, 0.12, nullptr);
  auto* c2 = new Cell(k, i, j, 0.53, nullptr);
  auto* c3 = new Cell(k, i, j, 0.71, nullptr);
  auto* c4 = new Cell(k, i, j, 0.88, nullptr);
  expected_tbl(k, i, j).push_back(c1);
  expected_tbl(k, i, j).push_back(c2);
  expected_tbl(k, i, j).push_back(c3);
  expected_tbl(k, i, j).push_back(c4);
  auto* c5 = new Cell(kp, ip, jp, 0.89, c1);
  auto* c6 = new Cell(kp, ip, jp, 1.3, c2);
  auto* c7 = new Cell(kp, ip, jp, 1.48, c3);
  auto* c8 = new Cell(kp, ip, jp, 1.65, c4);
  expected_tbl(kp, ip, jp).push_back(c5);
  expected_tbl(kp, ip, jp).push_back(c6);
  expected_tbl(kp, ip, jp).push_back(c7);
  expected_tbl(kp, ip, jp).push_back(c8);

  {
    nbest::NBestCells& cells = tbl(k, i, j);
    std::sort(cells.begin(), cells.end(), nbest::CellCompWithScoreReversed);
    nbest::NBestCells& ecells = expected_tbl(k, i, j);
    std::sort(ecells.begin(), ecells.end(), nbest::CellCompWithScoreReversed);
    EXPECT_EQ(cells.size(), ecells.size());
    EXPECT_TRUE(nbest::Equiv(expected_tbl(k, i, j), tbl(k, i, j)));
    // for_each(cells.begin(), cells.end(), nbest::ObjectDeleter());
    // for_each(ecells.begin(), ecells.end(), nbest::ObjectDeleter());
  }

  {
    nbest::NBestCells& cells = tbl(kp, ip, jp);
    std::sort(cells.begin(), cells.end(), nbest::CellCompWithScoreReversed);
    nbest::NBestCells& ecells = expected_tbl(kp, ip, jp);
    std::sort(ecells.begin(), ecells.end(), nbest::CellCompWithScoreReversed);
    EXPECT_EQ(cells.size(), ecells.size());
    EXPECT_TRUE(nbest::Equiv(expected_tbl(kp, ip, jp), tbl(kp, ip, jp)));
    // for_each(cells.begin(), cells.end(), nbest::ObjectDeleter());
    // for_each(ecells.begin(), ecells.end(), nbest::ObjectDeleter());
  }
}

// This test..
TEST_F(RecurrencesNbestTest, CellInserterBifurcation) {
  const int rna_size = 40;
  int k, i, j;
  std::tie(k, i, j) = std::make_tuple(1, 2, 3);
  int kp, ip, jp;
  std::tie(kp, ip, jp) = std::make_tuple(9, 22, 33);
  int kpp, ipp, jpp;
  std::tie(kpp, ipp, jpp) = std::make_tuple(4, 10, 5);
  const vector<double> v = {0.131, 0.177, 0.213, 0.311,
                            0.411, 0.194, 0.712, 0.121};

  nbest::DPTable tbl(NUM_STATES, rna_size, rna_size);
  auto* c1 = new Cell(k, i, j, v[0], nullptr);
  auto* c2 = new Cell(k, i, j, v[1], nullptr);
  auto* c3 = new Cell(k, i, j, v[2], nullptr);
  auto* c4 = new Cell(k, i, j, v[3], nullptr);
  auto* c5 = new Cell(kp, ip, jp, v[4], nullptr);
  auto* c6 = new Cell(kp, ip, jp, v[5], nullptr);
  auto* c7 = new Cell(kp, ip, jp, v[6], nullptr);
  auto* c8 = new Cell(kp, ip, jp, v[7], nullptr);
  tbl(k, i, j).push_back(c1);
  tbl(k, i, j).push_back(c2);
  tbl(k, i, j).push_back(c3);
  tbl(k, i, j).push_back(c4);
  tbl(kp, ip, jp).push_back(c5);
  tbl(kp, ip, jp).push_back(c6);
  tbl(kp, ip, jp).push_back(c7);
  tbl(kp, ip, jp).push_back(c8);
  // Insert new Cells
  const double s = 0.913;
  nbest::CreateCellsFromTwoCoordinates(&tbl, kpp, ipp, jpp, k, i, j, kp, ip, jp,
                                       s);

  nbest::DPTable exptbl(NUM_STATES, rna_size, rna_size);
  {
    nbest::NBestCells& t = exptbl(k, i, j);
    t.push_back(new Cell(k, i, j, v[0], nullptr));
    t.push_back(new Cell(k, i, j, v[1], nullptr));
    t.push_back(new Cell(k, i, j, v[2], nullptr));
    t.push_back(new Cell(k, i, j, v[3], nullptr));
  }
  {
    nbest::NBestCells& t = exptbl(kp, ip, jp);
    t.push_back(new Cell(kp, ip, jp, v[4], nullptr));
    t.push_back(new Cell(kp, ip, jp, v[5], nullptr));
    t.push_back(new Cell(kp, ip, jp, v[6], nullptr));
    t.push_back(new Cell(kp, ip, jp, v[7], nullptr));
  }
  {
    nbest::NBestCells& t = exptbl(kpp, ipp, jpp);
    t.push_back(new Cell(kpp, ipp, jpp, v[0] + v[4] + s, c1, c5));
    t.push_back(new Cell(kpp, ipp, jpp, v[0] + v[5] + s, c1, c6));
    t.push_back(new Cell(kpp, ipp, jpp, v[0] + v[6] + s, c1, c7));
    t.push_back(new Cell(kpp, ipp, jpp, v[0] + v[7] + s, c1, c8));
    t.push_back(new Cell(kpp, ipp, jpp, v[1] + v[4] + s, c2, c5));
    t.push_back(new Cell(kpp, ipp, jpp, v[1] + v[5] + s, c2, c6));
    t.push_back(new Cell(kpp, ipp, jpp, v[1] + v[6] + s, c2, c7));
    t.push_back(new Cell(kpp, ipp, jpp, v[1] + v[7] + s, c2, c8));
    t.push_back(new Cell(kpp, ipp, jpp, v[2] + v[4] + s, c3, c5));
    t.push_back(new Cell(kpp, ipp, jpp, v[2] + v[5] + s, c3, c6));
    t.push_back(new Cell(kpp, ipp, jpp, v[2] + v[6] + s, c3, c7));
    t.push_back(new Cell(kpp, ipp, jpp, v[2] + v[7] + s, c3, c8));
    t.push_back(new Cell(kpp, ipp, jpp, v[3] + v[4] + s, c4, c5));
    t.push_back(new Cell(kpp, ipp, jpp, v[3] + v[5] + s, c4, c6));
    t.push_back(new Cell(kpp, ipp, jpp, v[3] + v[6] + s, c4, c7));
    t.push_back(new Cell(kpp, ipp, jpp, v[3] + v[7] + s, c4, c8));
  }

  {
    nbest::NBestCells& cells = tbl(k, i, j);
    std::sort(cells.begin(), cells.end(), nbest::CellCompWithScoreReversed);
    nbest::NBestCells& ecells = exptbl(k, i, j);
    std::sort(ecells.begin(), ecells.end(), nbest::CellCompWithScoreReversed);
    EXPECT_EQ(cells.size(), ecells.size());
    EXPECT_TRUE(nbest::Equiv(ecells, cells));
  }
  {
    nbest::NBestCells& cells = tbl(kp, ip, jp);
    std::sort(cells.begin(), cells.end(), nbest::CellCompWithScoreReversed);
    nbest::NBestCells& ecells = exptbl(kp, ip, jp);
    std::sort(ecells.begin(), ecells.end(), nbest::CellCompWithScoreReversed);
    EXPECT_EQ(cells.size(), ecells.size());
    EXPECT_TRUE(nbest::Equiv(ecells, cells));
  }
  {
    nbest::NBestCells& cells = tbl(kpp, ipp, jpp);
    std::sort(cells.begin(), cells.end(), nbest::CellCompWithScoreReversed);
    nbest::NBestCells& ecells = exptbl(kpp, ipp, jpp);
    std::sort(ecells.begin(), ecells.end(), nbest::CellCompWithScoreReversed);
    EXPECT_EQ(cells.size(), ecells.size());
    EXPECT_TRUE(nbest::Equiv(ecells, cells));
  }
}

//
TEST_F(RecurrencesNbestTest, FillNBestTable) {
  Structure rna("(((..((...))....)))", "ACCATCTCCCGAACCAGCT", {});
  const int rna_size = rna.Size();
  FeatureVec model;
  const double s = 12.34;
  std::default_random_engine gen(s);
  for (int i = 0; i < Features_RNAFeatures_FV_SIZE; i++) {
    model.insert(make_pair(i, gen()));
  }
  InitializeFeatures(&model);
  // Create and fill dynamic programming tables.
  nbest::DPTable tbl(NUM_STATES, rna_size, rna_size);
  nbest::FillTable(rna, LOSS_ENABLED, 1.0, model, &tbl);
}

TEST_F(RecurrencesNbestTest, DecodeNBest) {
  Structure rna("(((..((...))....)))", "ACCATCTCCCGAACCAGCT", {});
  const idx_t rna_size = rna.Size();
  FeatureVec model;
  const double s = 12.34;
  std::default_random_engine gen(s);
  for (int i = 0; i < Features_RNAFeatures_FV_SIZE; i++) {
    model.insert(make_pair(i, static_cast<double>(gen()) / gen.max()));
  }
  InitializeFeatures(&model);
  const DECODING_MODE vmode = LOSS_DISABLED;
  const score_t loss_factor = 1.0;

  // Create one-best dp.
  DPTable* scores = new DPTable(NUM_STATES, rna_size, rna_size, LOG_ZERO);
  DPTable* traceback = new DPTable(NUM_STATES, rna_size, rna_size, LOG_ZERO);
  FillTables(rna, vmode, loss_factor, model, scores, traceback);
  vector<idx_t> result;
  const double best_path_score = scores->Get(DO_OUTER, 0, rna_size - 1);
  LOG(INFO) << "One-best: " << best_path_score;
  PerformTraceback(*traceback, rna_size, &result);
  string current;
  PairingsToBrackets(result, &current);
  LOG(INFO) << current;
  delete scores;
  delete traceback;

  // Create and fill dynamic programming tables.
  nbest::DPTable tbl(NUM_STATES, rna_size, rna_size);
  nbest::FillTable(rna, vmode, loss_factor, model, &tbl);

  const int n = 10;
  vector<vector<idx_t>> results;
  vector<double> nbest_scores;
  LOG(INFO) << "Scores are: ";
  for (int i = 0; i < n; i++) {
    vector<idx_t> res;
    nbest_scores.push_back(nbest::PerformTraceback(i, tbl, rna_size, &res));
    LOG(INFO) << nbest_scores.back();
    results.push_back(res);
    PairingsToBrackets(res, &current);
    LOG(INFO) << current;
  }
  EXPECT_LT(fabs(best_path_score - nbest_scores[0]), epsilon);
}

TEST_F(RecurrencesNbestTest, DecodeNBestLossEnabled) {
  Structure rna("(((..((...))....)))", "ACCATCTCCCGAACCAGCT", {});
  const idx_t rna_size = rna.Size();
  FeatureVec model;
  const double s = 12.34;
  std::default_random_engine gen(s);
  for (int i = 0; i < Features_RNAFeatures_FV_SIZE; i++) {
    model.insert(make_pair(i, static_cast<double>(gen()) / gen.max()));
  }
  InitializeFeatures(&model);
  const DECODING_MODE vmode = LOSS_ENABLED;
  const score_t loss_factor = 1.0;

  // Create one-best dp.
  DPTable* scores = new DPTable(NUM_STATES, rna_size, rna_size, LOG_ZERO);
  DPTable* traceback = new DPTable(NUM_STATES, rna_size, rna_size, LOG_ZERO);
  FillTables(rna, vmode, loss_factor, model, scores, traceback);
  vector<idx_t> result;
  const double best_path_score = scores->Get(DO_OUTER, 0, rna_size - 1);
  LOG(INFO) << "One-best: " << best_path_score;
  PerformTraceback(*traceback, rna_size, &result);
  string current;
  PairingsToBrackets(result, &current);
  LOG(INFO) << current;
  delete scores;
  delete traceback;

  // Create and fill dynamic programming tables.
  nbest::DPTable tbl(NUM_STATES, rna_size, rna_size);
  nbest::FillTable(rna, vmode, loss_factor, model, &tbl);

  const int n = 10;
  vector<vector<idx_t>> results;
  vector<double> nbest_scores;
  LOG(INFO) << "Scores are: ";
  for (int i = 0; i < n; i++) {
    vector<idx_t> res;
    nbest_scores.push_back(nbest::PerformTraceback(i, tbl, rna_size, &res));
    LOG(INFO) << nbest_scores.back();
    results.push_back(res);
    PairingsToBrackets(res, &current);
    LOG(INFO) << current;
  }
  EXPECT_LT(fabs(best_path_score - nbest_scores[0]), epsilon);
}

}  // namespace
}  // namespace qrsp
}  // namespace fdb
