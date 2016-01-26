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

#include <algorithm>
#include <cmath>
#include <string>

#include <gtest/gtest.h>
#include "learning-utils.h"
#include "model.h"
#include "third_party/googleflags/include/gflags/gflags.h"
#include "utils.h"

DECLARE_bool(enable_quality_features);

namespace qrisp {
namespace {

using fdb::learning::CalculateSparseScalarProduct;

constexpr const float kEpsilon = 1e-6;

class RecurrencesTest : public testing::Test {
 protected:
  RecurrencesTest() {}
  ~RecurrencesTest() {}
  // virtual void SetUp() will be called before each test is run. You
  // should define it if you need to initialize the varaibles.
  // Otherwise, this can be skipped.
  virtual void SetUp() {}

  virtual void TearDown() {}
};

TEST_F(RecurrencesTest, InitTables) {
  uint32_t table_size = 10;
  // Create and fill dynamic programming tables.
  DPTable score(NUM_STATES, table_size, table_size, LOG_ZERO);
  DPTable traceback(NUM_STATES, table_size, table_size, LOG_ZERO);
  for (idx_t i = 0; i < NUM_STATES; i++) {
    for (idx_t j = 0; j < table_size; j++) {
      for (idx_t k = 0; k < table_size; k++) {
        EXPECT_NEAR(score(i, j, k), LOG_ZERO, kEpsilon);
        EXPECT_NEAR(traceback(i, j, k), LOG_ZERO, kEpsilon);
      }
    }
  }
}

TEST_F(RecurrencesTest, UpdateTables) {
  idx_t table_size = 4;
  // Create and fill dynamic programming tables.
  DPTable score(NUM_STATES, table_size, table_size, LOG_ZERO);
  DPTable traceback(NUM_STATES, table_size, table_size, LOG_ZERO);
  for (idx_t i = 0; i < NUM_STATES; i++) {
    for (idx_t j = 0; j < table_size; j++) {
      for (idx_t k = 0; k < table_size; k++) {
        EXPECT_NEAR(score(i, j, k), LOG_ZERO, kEpsilon);
        EXPECT_NEAR(traceback(i, j, k), LOG_ZERO, kEpsilon);
      }
    }
  }
  for (idx_t i = 0; i < NUM_STATES; i++) {
    for (idx_t j = 0; j < table_size; j++) {
      for (idx_t k = 0; k < table_size; k++) {
        double new_score = (i + j) % 2 == 0 ? k * 1.0 : k * -1.0;
        double new_trace = (i + j) % 2 == 0 ? 55.0 : 2.0;
        UpdateMax(&score, &traceback, i, j, k, new_score, new_trace);
      }
    }
  }
  for (idx_t i = 0; i < NUM_STATES; i++) {
    for (idx_t j = 0; j < table_size; j++) {
      for (idx_t k = 0; k < table_size; k++) {
        EXPECT_NEAR(score(i, j, k), (i + j) % 2 == 0 ? k * 1.0 : k * -1.0,
                    kEpsilon);
        EXPECT_NEAR(traceback(i, j, k), (i + j) % 2 == 0 ? 55.0 : 2.0,
                    kEpsilon);
      }
    }
  }
  for (idx_t i = 0; i < NUM_STATES; i++) {
    for (idx_t j = 0; j < table_size; j++) {
      for (idx_t k = 0; k < table_size; k++) {
        double new_score = (i + j) % 2 == 0 ? k * 10.0 : k * -0.1;
        double new_trace = (i + j) % 2 == 0 ? 33.0 : 4.0;
        UpdateMax(&score, &traceback, i, j, k, new_score, new_trace);
      }
    }
  }
  for (idx_t i = 0; i < NUM_STATES; i++) {
    for (idx_t j = 0; j < table_size; j++) {
      for (idx_t k = 0; k < table_size; k++) {
        if (k > 0) {
          EXPECT_NEAR(score(i, j, k), (i + j) % 2 == 0 ? k * 10.0 : k * -0.1,
                      kEpsilon);
          EXPECT_NEAR(traceback(i, j, k), (i + j) % 2 == 0 ? 33.0 : 4.0,
                      kEpsilon);
        } else {
          EXPECT_NEAR(score(i, j, k), (i + j) % 2 == 0 ? k * 1.0 : k * -1.0,
                      kEpsilon);
          EXPECT_NEAR(traceback(i, j, k), (i + j) % 2 == 0 ? 55.0 : 2.0,
                      kEpsilon);
        }
      }
    }
  }
}

TEST_F(RecurrencesTest, UpdateTablesViaFunctor) {
  uint32_t table_size = 4;
  // Create and fill dynamic programming tables.
  DPTable score(NUM_STATES, table_size, table_size, LOG_ZERO);
  DPTable traceback(NUM_STATES, table_size, table_size, LOG_ZERO);
  for (idx_t i = 0; i < NUM_STATES; i++) {
    for (idx_t j = 0; j < table_size; j++) {
      for (idx_t k = 0; k < table_size; k++) {
        EXPECT_NEAR(score(i, j, k), LOG_ZERO, kEpsilon);
        EXPECT_NEAR(traceback(i, j, k), LOG_ZERO, kEpsilon);
      }
    }
  }

  auto UPDATE = bind(UpdateMax, &score, &traceback, std::placeholders::_1,
                     std::placeholders::_2, std::placeholders::_3,
                     std::placeholders::_4, std::placeholders::_5);

  for (idx_t i = 0; i < NUM_STATES; i++) {
    for (idx_t j = 0; j < table_size; j++) {
      for (idx_t k = 0; k < table_size; k++) {
        double new_score = (i + j) % 2 == 0 ? k * 1.0 : k * -1.0;
        double new_trace = (i + j) % 2 == 0 ? 55.0 : 2.0;
        UPDATE(i, j, k, new_score, new_trace);
      }
    }
  }
  for (idx_t i = 0; i < NUM_STATES; i++) {
    for (idx_t j = 0; j < table_size; j++) {
      for (idx_t k = 0; k < table_size; k++) {
        EXPECT_NEAR(score(i, j, k), (i + j) % 2 == 0 ? k * 1.0 : k * -1.0,
                    kEpsilon);
        EXPECT_NEAR(traceback(i, j, k), (i + j) % 2 == 0 ? 55.0 : 2.0,
                    kEpsilon);
      }
    }
  }
}

TEST_F(RecurrencesTest, LossFunctor) {}

TEST_F(RecurrencesTest, EncodeDecodeTraceback) {
  const idx_t rna_size = 19;
  {
    idx_t i = 3;
    idx_t j = 15;
    score_t enc_tb = EncodeTraceback(i, j, rna_size);
    idx_t ip;
    idx_t jp;
    std::tie(ip, jp) = DecodeTraceback(enc_tb, rna_size);
    EXPECT_EQ(i, ip);
    EXPECT_EQ(j, jp);
  }

  {
    score_t enc_tb = EncodeTraceback(0, 0, rna_size);
    idx_t ip;
    idx_t jp;
    std::tie(ip, jp) = DecodeTraceback(enc_tb, rna_size);
    EXPECT_EQ(0, ip);
    EXPECT_EQ(0, jp);
  }

  {
    score_t enc_tb = EncodeTraceback(1, 0, rna_size);
    idx_t ip;
    idx_t jp;
    std::tie(ip, jp) = DecodeTraceback(enc_tb, rna_size);
    EXPECT_EQ(1, ip);
    EXPECT_EQ(0, jp);
  }

  const idx_t rna_size_2 = 12;
  {
    score_t enc_tb = EncodeTraceback(2, 0, rna_size_2);
    idx_t ip;
    idx_t jp;
    std::tie(ip, jp) = DecodeTraceback(enc_tb, rna_size_2);
    EXPECT_EQ(2, ip);
    EXPECT_EQ(0, jp);
  }

  {
    score_t enc_tb = EncodeTraceback(0, 0, rna_size_2);
    idx_t ip;
    idx_t jp;
    std::tie(ip, jp) = DecodeTraceback(enc_tb, rna_size_2);
    EXPECT_EQ(0, ip);
    EXPECT_EQ(0, jp);
  }

  {
    score_t enc_tb = EncodeTraceback(1, 0, rna_size_2);
    idx_t ip;
    idx_t jp;
    std::tie(ip, jp) = DecodeTraceback(enc_tb, rna_size_2);
    EXPECT_EQ(1, ip);
    EXPECT_EQ(0, jp);
  }

  {
    score_t enc_tb = EncodeTraceback(2, 0, rna_size);
    idx_t ip;
    idx_t jp;
    std::tie(ip, jp) = DecodeTraceback(enc_tb, rna_size);
    EXPECT_EQ(2, ip);
    EXPECT_EQ(0, jp);
  }

  {
    score_t enc_tb = EncodeTraceback(3, 11, rna_size);
    idx_t ip;
    idx_t jp;
    std::tie(ip, jp) = DecodeTraceback(enc_tb, rna_size);
    EXPECT_EQ(3, ip);
    EXPECT_EQ(11, jp);
  }
}

TEST_F(RecurrencesTest, FillTables) {
  Structure rna(".......", "ACTCACT", {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
  const idx_t rna_size = rna.Size();
  // Create and fill dynamic programming tables.
  DPTable* scores = new DPTable(NUM_STATES, rna_size, rna_size, LOG_ZERO);
  DPTable* traceback = new DPTable(NUM_STATES, rna_size, rna_size, LOG_ZERO);
  DECODING_MODE vmode = LOSS_DISABLED;
  const score_t loss_factor = 1.0;
  FeatureVec model;
  CalculateFeatures(rna, &model);
  FillTables(rna, vmode, loss_factor, model, scores, traceback);
  delete scores;
  delete traceback;
}

TEST_F(RecurrencesTest, PerformTraceback) {
  vector<score_t> q(19, 0.0);
  Structure rna("...(((....)))......", "ACCATCTCCCGAACCAGCT", q);
  const idx_t rna_size = rna.Size();
  DPTable* scores = new DPTable(NUM_STATES, rna_size, rna_size, LOG_ZERO);
  DPTable* traceback = new DPTable(NUM_STATES, rna_size, rna_size, LOG_ZERO);
  const DECODING_MODE vmode = LOSS_DISABLED;
  const score_t loss_factor = 1.0;
  FeatureVec model;
  CalculateFeatures(rna, &model);
  FillTables(rna, vmode, loss_factor, model, scores, traceback);
  vector<idx_t> result;
  PerformTraceback(*traceback, rna_size, &result);
  delete scores;
  delete traceback;
  // EXPECT_NEAR(196.0, scores.Get(DO_OUTER, 0, rna_size - 1), kEpsilon);
}

TEST_F(RecurrencesTest, DecodeBestPath) {
  vector<score_t> q(19, 0.0);
  Structure rna("...(((....)))......", "ACCATCTCCCGAACCAGCT", q);
  const idx_t rna_size = rna.Size();
  DPTable* scores = new DPTable(NUM_STATES, rna_size, rna_size, LOG_ZERO);
  DPTable* traceback = new DPTable(NUM_STATES, rna_size, rna_size, LOG_ZERO);
  FeatureVec model;
  CalculateFeatures(rna, &model);
  const DECODING_MODE vmode = LOSS_DISABLED;
  const score_t loss_factor = 1.0;
  FillTables(rna, vmode, loss_factor, model, scores, traceback);
  vector<idx_t> result;
  PerformTraceback(*traceback, rna_size, &result);
  string brackets;
  PairingsToBrackets(result, &brackets);
  Structure dec_rna(brackets, "ACCATCTCCCGAACCAGCT", q);
  FeatureVec features;
  CalculateFeatures(dec_rna, &features);
  // cout << "score best path: " << scores->Get(DO_OUTER, 0, rna_size - 1) <<
  // endl;
  // cout << "score scalar product: " << CalculateSparseScalarProduct(model,
  // features) << endl;
  delete scores;
  delete traceback;
}

//  for (int s = 0; s < NUM_STATES; s++) {
//    for (int i = 0; i < traceback->GetSize2(); i++) {
//      for (int j = 0; j < traceback->GetSize3(); j++) {
//        if (traceback->Get(s, i, j) >= 0) {
//          int p, k;
//          std::tie(p, k) = DecodeTraceback(traceback->Get(s, i, j), rna_size);
//          ASSERT_GE(p, 0);
//          ASSERT_GE(k, 0);
//          ASSERT_LE(p, rna_size);
//          ASSERT_LE(k, rna_size);
//        }
//      }
//    }
//  }

TEST_F(RecurrencesTest, DecodeBestPathRandShort) {
  Structure rna(".......", "ACCCTTA", {});
  const idx_t rna_size = rna.Size();
  DPTable* scores = new DPTable(NUM_STATES, rna_size, rna_size, LOG_ZERO);
  DPTable* traceback = new DPTable(NUM_STATES, rna_size, rna_size, LOG_ZERO);
  FeatureVec model;
  int ctr = 1;
  auto GetValue = [&ctr]() { return drand48(); };
  for (int i = 0; i < 906; i++) {
    model.insert(make_pair(i, GetValue()));
  }
  const DECODING_MODE vmode = LOSS_DISABLED;
  const score_t loss_factor = 1.0;
  FeatureVec model_transformed = model;
  InitializeFeatures(&model_transformed);
  FillTables(rna, vmode, loss_factor, model_transformed, scores, traceback);
  vector<idx_t> result;
  PerformTraceback(*traceback, rna_size, &result);
  vector<idx_t> index(result.size(), 0);
  int count = 0;
  for_each(index.begin(), index.end(),
           [&count](idx_t& i) { return i = count++; });
  // for (int i = 1; i < index.size(); i++) {
  //  cout << setw(3) << index[i] << " ";
  //}
  // cout << endl;
  // for (int i = 1; i < index.size(); i++) {
  //  cout << setw(3) << result[i] << " ";
  //}
  // cout << endl;
  string brackets;
  PairingsToBrackets(result, &brackets);
  Structure dec_rna(brackets, "ACCCTTA", {});
  FeatureVec features;
  CalculateFeatures(dec_rna, &features);
  EXPECT_EQ(scores->Get(DO_OUTER, 0, rna_size - 1),
            CalculateSparseScalarProduct(model, features));
  delete scores;
  delete traceback;
}

TEST_F(RecurrencesTest, DecodeBestPathRand) {
  // Structure rna("...(((....)))......",
  //                 "ACCATCTCCCGAACCAGCT");
  Structure rna("............", "ACCTTTTTTATT", {});
  const idx_t rna_size = rna.Size();
  DPTable* scores = new DPTable(NUM_STATES, rna_size, rna_size, LOG_ZERO);
  DPTable* traceback = new DPTable(NUM_STATES, rna_size, rna_size, LOG_ZERO);
  FeatureVec model;
  int ctr = 1;
  auto GetValue = [&ctr]() { return drand48(); };
  for (int i = 0; i < 906; i++) {
    model.insert(make_pair(i, GetValue()));
  }
  const DECODING_MODE vmode = LOSS_DISABLED;
  const score_t loss_factor = 1.0;
  FeatureVec model_transformed = model;
  InitializeFeatures(&model_transformed);
  FillTables(rna, vmode, loss_factor, model_transformed, scores, traceback);
  vector<idx_t> result;
  PerformTraceback(*traceback, rna_size, &result);
  vector<idx_t> index(result.size(), 0);
  int count = 0;
  for_each(index.begin(), index.end(),
           [&count](idx_t& i) { return i = count++; });
  // for (int i = 1; i < index.size(); i++) {
  //  cout << setw(3) << index[i] << " ";
  //}
  // cout << endl;
  // for (int i = 1; i < index.size(); i++) {
  //  cout << setw(3) << result[i] << " ";
  //}
  // cout << endl;
  string brackets;
  PairingsToBrackets(result, &brackets);
  // Structure dec_rna(brackets, "ACCATCTCCCGAACCAGCT");
  Structure dec_rna(brackets, "ACCTTTTTTATT", {});
  FeatureVec features;
  CalculateFeatures(dec_rna, &features);
  // cout << "Brackets: " << brackets << endl;
  // cout << "Score best path / scalar prod.: " << scores->Get(DO_OUTER, 0,
  // rna_size - 1)
  //     << " "
  //     << CalculateSparseScalarProduct(model, features)
  //     << endl;
  EXPECT_NEAR(scores->Get(DO_OUTER, 0, rna_size - 1),
              CalculateSparseScalarProduct(model, features), kEpsilon);
  delete scores;
  delete traceback;
}

TEST_F(RecurrencesTest, DecodeBestPathRandWithLoss) {
  // Structure rna("...(((....)))......",
  //                 "ACCATCTCCCGAACCAGCT");
  Structure rna("............", "ACCTTTTTTATT", {});
  const idx_t rna_size = rna.Size();
  DPTable* scores = new DPTable(NUM_STATES, rna_size, rna_size, LOG_ZERO);
  DPTable* traceback = new DPTable(NUM_STATES, rna_size, rna_size, LOG_ZERO);
  FeatureVec model;
  int ctr = 1;
  auto GetValue = [&ctr]() { return drand48(); };
  for (int i = 0; i < 906; i++) {
    model.insert(make_pair(i, GetValue()));
  }
  const DECODING_MODE vmode = LOSS_ENABLED;
  const score_t loss_factor = 1.0;
  FeatureVec model_transformed = model;
  InitializeFeatures(&model_transformed);
  FillTables(rna, vmode, loss_factor, model_transformed, scores, traceback);
  vector<idx_t> result;
  PerformTraceback(*traceback, rna_size, &result);
  vector<idx_t> index(result.size(), 0);
  int count = 0;
  for_each(index.begin(), index.end(),
           [&count](idx_t& i) { return i = count++; });
  for (idx_t i = 1; i < index.size(); i++) {
    cout << setw(3) << index[i] << " ";
  }
  cout << endl;
  for (idx_t i = 1; i < index.size(); i++) {
    cout << setw(3) << result[i] << " ";
  }
  cout << endl;
  string brackets;
  PairingsToBrackets(result, &brackets);
  // Structure dec_rna(brackets, "ACCATCTCCCGAACCAGCT");
  Structure dec_rna(brackets, "ACCTTTTTTATT", {});
  score_t loss = HammingLoss(rna, dec_rna);
  FeatureVec features;
  CalculateFeatures(dec_rna, &features);
  // cout << "Brackets: " << brackets << endl;
  // cout << "Score best path / scalar prod.: " << scores->Get(DO_OUTER, 0,
  // rna_size - 1)
  //     << " "
  //     << CalculateSparseScalarProduct(model, features) + loss
  //     << endl;
  EXPECT_NEAR(scores->Get(DO_OUTER, 0, rna_size - 1),
              CalculateSparseScalarProduct(model, features) + loss, kEpsilon);
  delete scores;
  delete traceback;
}

TEST_F(RecurrencesTest, DecodeBestPathRandLarge) {
  Structure rna("...(((....)))......", "ACCATCTCCCGAACCAGCT", {});
  const idx_t rna_size = rna.Size();
  DPTable* scores = new DPTable(NUM_STATES, rna_size, rna_size, LOG_ZERO);
  DPTable* traceback = new DPTable(NUM_STATES, rna_size, rna_size, LOG_ZERO);
  FeatureVec model;
  int ctr = 1;
  auto GetValue = [&ctr]() { return drand48(); };
  for (idx_t i = 0; i < 906; i++) {
    model.insert(make_pair(i, GetValue()));
  }
  const DECODING_MODE vmode = LOSS_DISABLED;
  const score_t loss_factor = 1.0;
  FeatureVec model_transformed = model;
  InitializeFeatures(&model_transformed);
  FillTables(rna, vmode, loss_factor, model_transformed, scores, traceback);
  vector<idx_t> result;
  PerformTraceback(*traceback, rna_size, &result);
  vector<idx_t> index(result.size(), 0);
  int count = 0;
  for_each(index.begin(), index.end(),
           [&count](idx_t& i) { return i = count++; });
  // for (int i = 1; i < index.size(); i++) {
  //  cout << setw(3) << index[i] << " ";
  //}
  // cout << endl;
  // for (int i = 1; i < index.size(); i++) {
  //  cout << setw(3) << result[i] << " ";
  //}
  // cout << endl;
  string brackets;
  PairingsToBrackets(result, &brackets);
  cout << brackets << endl;
  Structure dec_rna(brackets, "ACCATCTCCCGAACCAGCT", {});
  FeatureVec features;
  CalculateFeatures(dec_rna, &features);
  // cout << "Brackets: " << brackets << endl;
  // cout << "Score best path / scalar prod.: " << scores->Get(DO_OUTER, 0,
  // rna_size - 1)
  //     << " "
  //     << CalculateSparseScalarProduct(model, features)
  //     << endl;
  EXPECT_NEAR(scores->Get(DO_OUTER, 0, rna_size - 1),
              CalculateSparseScalarProduct(model, features), kEpsilon);
  delete scores;
  delete traceback;
}

TEST_F(RecurrencesTest, DecodeBestPathRandLargeLoss) {
  Structure rna("...(((....)))......", "ACCATCTCCCGAACCAGCT", {});
  const idx_t rna_size = rna.Size();
  DPTable* scores = new DPTable(NUM_STATES, rna_size, rna_size, LOG_ZERO);
  DPTable* traceback = new DPTable(NUM_STATES, rna_size, rna_size, LOG_ZERO);
  FeatureVec model;
  int ctr = 1;
  auto GetValue = [&ctr]() { return drand48(); };
  for (idx_t i = 0; i < 906; i++) {
    model.insert(make_pair(i, GetValue()));
  }
  const DECODING_MODE vmode = LOSS_ENABLED;
  const score_t loss_factor = 1.0;
  FeatureVec model_transformed = model;
  InitializeFeatures(&model_transformed);
  FillTables(rna, vmode, loss_factor, model_transformed, scores, traceback);
  vector<idx_t> result;
  PerformTraceback(*traceback, rna_size, &result);
  vector<idx_t> index(result.size(), 0);
  int count = 0;
  for_each(index.begin(), index.end(),
           [&count](idx_t& i) { return i = count++; });
  // for (int i = 1; i < index.size(); i++) {
  //  cout << setw(3) << index[i] << " ";
  //}
  // cout << endl;
  // for (int i = 1; i < index.size(); i++) {
  //  cout << setw(3) << result[i] << " ";
  //}
  // cout << endl;
  string brackets;
  PairingsToBrackets(result, &brackets);
  Structure dec_rna(brackets, "ACCATCTCCCGAACCAGCT", {});
  score_t loss = HammingLoss(rna, dec_rna);
  FeatureVec features;
  CalculateFeatures(dec_rna, &features);
  // cout << "Brackets: " << brackets << endl;
  // cout << "Score best path / scalar prod.: " << scores->Get(DO_OUTER, 0,
  // rna_size - 1)
  //     << " "
  //     << CalculateSparseScalarProduct(model, features) + loss
  //     << endl;
  EXPECT_NEAR(scores->Get(DO_OUTER, 0, rna_size - 1),
              CalculateSparseScalarProduct(model, features) + loss, kEpsilon);
  delete scores;
  delete traceback;
}

}  // namespace
}  // namespace qrisp

int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
