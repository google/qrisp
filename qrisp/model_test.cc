#include "model.h"
#include "utils.h"

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <sstream>
#include <string>
#include "qrisp/proto/parameters.pb.h"
#include "structure.h"
#include "third_party/googleflags/include/gflags/gflags.h"

namespace qrisp {
namespace {

constexpr const float kEpsilon = 1e-6;

using std::make_tuple;

class ModelTest : public testing::Test {
 protected:
  ModelTest() {}
  ~ModelTest() {}

  virtual void SetUp() {}
  virtual void TearDown() {}
};

TEST_F(ModelTest, CalculateFullOffset) {
  EXPECT_EQ(0, CalculateFullOffset(make_tuple(0, 0, 0, 0)));
  EXPECT_EQ(85, CalculateFullOffset(make_tuple(1, 1, 1, 1)));
  EXPECT_EQ(170, CalculateFullOffset(make_tuple(2, 2, 2, 2)));
  EXPECT_EQ(255, CalculateFullOffset(make_tuple(3, 3, 3, 3)));

  EXPECT_EQ(1, CalculateFullOffset(make_tuple(0, 0, 0, 1)));
  EXPECT_EQ(64, CalculateFullOffset(make_tuple(1, 0, 0, 0)));
  EXPECT_EQ(2, CalculateFullOffset(make_tuple(0, 0, 0, 2)));
  EXPECT_EQ(128, CalculateFullOffset(make_tuple(2, 0, 0, 0)));
  EXPECT_EQ(3, CalculateFullOffset(make_tuple(0, 0, 0, 3)));
  EXPECT_EQ(192, CalculateFullOffset(make_tuple(3, 0, 0, 0)));
}

TEST_F(ModelTest, CalculateOffset) {
  EXPECT_EQ(0, CalculateOffset(make_tuple(0, 0, 0, 0)));
  EXPECT_EQ(72, CalculateOffset(make_tuple(1, 1, 1, 1)));
  EXPECT_EQ(117, CalculateOffset(make_tuple(2, 2, 2, 2)));
  EXPECT_EQ(135, CalculateOffset(make_tuple(3, 3, 3, 3)));

  EXPECT_EQ(1, CalculateOffset(make_tuple(1, 0, 0, 0)));
  EXPECT_EQ(1, CalculateOffset(make_tuple(0, 0, 0, 1)));
  EXPECT_EQ(2, CalculateOffset(make_tuple(2, 0, 0, 0)));
  EXPECT_EQ(2, CalculateOffset(make_tuple(0, 0, 0, 2)));
  EXPECT_EQ(3, CalculateOffset(make_tuple(3, 0, 0, 0)));
  EXPECT_EQ(3, CalculateOffset(make_tuple(0, 0, 0, 3)));

  EXPECT_EQ(127, CalculateOffset(make_tuple(3, 0, 1, 3)));
}

TEST_F(ModelTest, CalculateLeftOffset) {
  EXPECT_EQ(0, CalculateLeftOffset(make_tuple(0, 0, 0, IDX_NOT_SET)));
  EXPECT_EQ(1, CalculateLeftOffset(make_tuple(0, 0, 1, IDX_NOT_SET)));

  EXPECT_EQ(20, CalculateLeftOffset(make_tuple(1, 1, 0, IDX_NOT_SET)));
  EXPECT_EQ(21, CalculateLeftOffset(make_tuple(1, 1, 1, IDX_NOT_SET)));
  EXPECT_EQ(22, CalculateLeftOffset(make_tuple(1, 1, 2, IDX_NOT_SET)));

  EXPECT_EQ(41, CalculateLeftOffset(make_tuple(2, 2, 1, IDX_NOT_SET)));
  EXPECT_EQ(42, CalculateLeftOffset(make_tuple(2, 2, 2, IDX_NOT_SET)));
  EXPECT_EQ(43, CalculateLeftOffset(make_tuple(2, 2, 3, IDX_NOT_SET)));

  EXPECT_EQ(62, CalculateLeftOffset(make_tuple(3, 3, 2, IDX_NOT_SET)));
  EXPECT_EQ(63, CalculateLeftOffset(make_tuple(3, 3, 3, IDX_NOT_SET)));
}

TEST_F(ModelTest, CalculateRightOffset) {
  EXPECT_EQ(0, CalculateRightOffset(make_tuple(0, 0, IDX_NOT_SET, 0)));
  EXPECT_EQ(1, CalculateRightOffset(make_tuple(0, 0, IDX_NOT_SET, 1)));

  EXPECT_EQ(20, CalculateRightOffset(make_tuple(1, 1, IDX_NOT_SET, 0)));
  EXPECT_EQ(21, CalculateRightOffset(make_tuple(1, 1, IDX_NOT_SET, 1)));
  EXPECT_EQ(22, CalculateRightOffset(make_tuple(1, 1, IDX_NOT_SET, 2)));

  EXPECT_EQ(41, CalculateRightOffset(make_tuple(2, 2, IDX_NOT_SET, 1)));
  EXPECT_EQ(42, CalculateRightOffset(make_tuple(2, 2, IDX_NOT_SET, 2)));
  EXPECT_EQ(43, CalculateRightOffset(make_tuple(2, 2, IDX_NOT_SET, 3)));

  EXPECT_EQ(62, CalculateRightOffset(make_tuple(3, 3, IDX_NOT_SET, 2)));
  EXPECT_EQ(63, CalculateRightOffset(make_tuple(3, 3, IDX_NOT_SET, 3)));
}

TEST_F(ModelTest, GetSummedParameters) {
  FeatureVec scores = {mp(0, 0.0), mp(1, 0.2),  mp(2, 0.3),
                       mp(3, 0.0), mp(4, -0.2), mp(5, 0.0)};
  IndexVec idx;
  // Check for no op.
  EXPECT_NEAR(GetSummedParameters(idx, scores), 0.0, kEpsilon);
  idx.push_back(mp(1, 0.2));
  idx.push_back(mp(2, 0.2));
  // Now check for some real feature summation.
  EXPECT_NEAR(GetSummedParameters(idx, scores), 0.10, kEpsilon);
}

TEST_F(ModelTest, UpdateFeatureWeights) {
  FeatureVec scores = {mp(0, 0.0), mp(1, 0.2),  mp(2, 0.3),
                       mp(3, 0.0), mp(4, -0.2), mp(5, 1.0)};
  IndexVec idx;
  const FeatureVec original_scores = scores;
  UpdateFeatureWeights(idx, &scores);
  // Check for no op.
  EXPECT_EQ(original_scores, scores);
  // Now check for some real updates.
  IndexVec idx2 = {mp(1, 0.2), mp(2, 0.2), mp(5, -0.2), mp(6, 0.77)};
  UpdateFeatureWeights(idx2, &scores);
  FeatureVec expected = {mp(0, 0.0),  mp(1, 0.4), mp(2, 0.5), mp(3, 0.0),
                         mp(4, -0.2), mp(5, 0.8), mp(6, 0.77)};
  EXPECT_EQ(expected, scores);
}

TEST_F(ModelTest, HairpinWeights) {
  // for (int i = 0; i < structures_.size(); i++) {
  //  const rna_t t = structures_[i];
  //  const Structure rna(get<0>(t), get<2>(t), get<4>(t));
  //  ScoreVec params;
  //  auto HairpinScorer = bind(GenericScorer, &HairpinWeights, rna,
  //                            placeholders::_1, params);
  //  auto HairpinCounter = bind(Counter, &HairpinScorer, placeholders::_1);

  //}
  // Tuple t(1, 2, 3, 4);
  // HairpinScorer_(t);
  // vector<score_t> vec(FV_SIZE, LOG_ZERO);
  // vector<Structure>::const_iterator it;
  // for (it = structures_.begin(); it != structures_.end(); ++it) {
  //  HairpinScorer* hairpin = new HairpinScorer(&rna, &vec, SCORE);
  //}
}

TEST_F(ModelTest, HelixBasePairWeights) {
  // enable_quality_features = false;
  vector<score_t> q(12, 0.0);
  Structure rna("............", "ACCACCTCCCGA", q);
  IndexVec indices;

  HelixBasePairWeights(make_tuple(1, 7, 1, 1), rna, &indices);
  IndexVec expected = {mp(Features::BASE_PAIR_AU, 1)};
  EXPECT_EQ(expected, indices);

  HelixBasePairWeights(make_tuple(4, 11, 1, 1), rna, &indices);
  expected = {mp(Features::BASE_PAIR_AU, 1), mp(Features::BASE_PAIR_AG, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  HelixBasePairWeights(make_tuple(4, 11, 1, 1), rna, &indices);
  expected = {mp(Features::BASE_PAIR_AG, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  HelixBasePairWeights(make_tuple(2, 6, 1, 1), rna, &indices);
  expected = {mp(Features::BASE_PAIR_CC, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  HelixBasePairWeights(make_tuple(1, 8, 1, 1), rna, &indices);
  expected = {mp(Features::BASE_PAIR_AC, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  HelixBasePairWeights(make_tuple(2, 12, 1, 1), rna, &indices);
  expected = {mp(Features::BASE_PAIR_AC, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  HelixBasePairWeights(make_tuple(2, 11, 1, 1), rna, &indices);
  expected = {mp(Features::BASE_PAIR_CG, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  HelixBasePairWeights(make_tuple(1, 4, 1, 1), rna, &indices);
  expected = {mp(Features::BASE_PAIR_AA, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();

  vector<score_t> q2(11, 0.0);
  Structure rna2("...........", "TGGTGGAGGGC", q2);
  HelixBasePairWeights(make_tuple(1, 7, 1, 1), rna2, &indices);
  expected = {mp(Features::BASE_PAIR_AU, 1)};
  EXPECT_EQ(expected, indices);

  HelixBasePairWeights(make_tuple(4, 11, 1, 1), rna2, &indices);
  expected = {mp(Features::BASE_PAIR_AU, 1), mp(Features::BASE_PAIR_CU, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  HelixBasePairWeights(make_tuple(4, 11, 1, 1), rna2, &indices);
  expected = {mp(Features::BASE_PAIR_CU, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  HelixBasePairWeights(make_tuple(2, 6, 1, 1), rna2, &indices);
  expected = {mp(Features::BASE_PAIR_GG, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  HelixBasePairWeights(make_tuple(1, 4, 1, 1), rna2, &indices);
  expected = {mp(Features::BASE_PAIR_UU, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  HelixBasePairWeights(make_tuple(1, 5, 1, 1), rna2, &indices);
  expected = {mp(Features::BASE_PAIR_GU, 1)};
  EXPECT_EQ(expected, indices);
}

TEST_F(ModelTest, HelixChangeWeights) {
  Structure rna("............", "ACCACCTCCCGA", {});
  IndexVec indices;

  HelixChangeWeights(make_tuple(1, 1, 1, 1), rna, &indices);
  IndexVec expected = {mp(Features::STACKING_1_SCORE, 1)};
  EXPECT_EQ(expected, indices);
  indices.clear();
  HelixChangeWeights(make_tuple(1, 1, 2, 1), rna, &indices);
  expected = {mp(Features::STACKING_2_SCORE, 1)};
  EXPECT_EQ(expected, indices);
  indices.clear();
  HelixChangeWeights(make_tuple(1, 1, 3, 1), rna, &indices);
  expected = {mp(Features::STACKING_3_SCORE, 1)};
  EXPECT_EQ(expected, indices);
  indices.clear();
  HelixChangeWeights(make_tuple(1, 1, 4, 1), rna, &indices);
  expected = {mp(Features::STACKING_4_SCORE, 1)};
  EXPECT_EQ(expected, indices);
  indices.clear();
  HelixChangeWeights(make_tuple(1, 1, 5, 1), rna, &indices);
  expected = {mp(Features::STACKING_5_SCORE, 1)};
  EXPECT_EQ(expected, indices);
  indices.clear();
  HelixChangeWeights(make_tuple(1, 1, 6, 1), rna, &indices);
  expected = {mp(Features::EXTENSION_SCORE, 1)};
  EXPECT_EQ(expected, indices);
  indices.clear();
  HelixChangeWeights(make_tuple(1, 1, 16, 1), rna, &indices);
  expected = {mp(Features::EXTENSION_SCORE, 1)};
  EXPECT_EQ(expected, indices);
}

TEST_F(ModelTest, HelixClosingWeights) {
  Structure rna("............", "ACCACCTCCCGA", {});
  IndexVec indices;

  HelixClosingWeights(make_tuple(1, 7, 1, 1), rna, &indices);
  IndexVec expected = {mp(Features::HELIX_CLOSING_AU, 1)};
  EXPECT_EQ(expected, indices);

  HelixClosingWeights(make_tuple(4, 11, 1, 1), rna, &indices);
  expected = {mp(Features::HELIX_CLOSING_AU, 1),
              mp(Features::HELIX_CLOSING_AG, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  HelixClosingWeights(make_tuple(4, 11, 1, 1), rna, &indices);
  expected = {mp(Features::HELIX_CLOSING_AG, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  HelixClosingWeights(make_tuple(2, 6, 1, 1), rna, &indices);
  expected = {mp(Features::HELIX_CLOSING_CC, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  HelixClosingWeights(make_tuple(1, 8, 1, 1), rna, &indices);
  expected = {mp(Features::HELIX_CLOSING_AC, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  HelixClosingWeights(make_tuple(2, 12, 1, 1), rna, &indices);
  expected = {mp(Features::HELIX_CLOSING_CA, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  HelixClosingWeights(make_tuple(2, 11, 1, 1), rna, &indices);
  expected = {mp(Features::HELIX_CLOSING_CG, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  HelixClosingWeights(make_tuple(1, 4, 1, 1), rna, &indices);
  expected = {mp(Features::HELIX_CLOSING_AA, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  Structure rna2("...........", "TGGTGGAGGGC", {});
  HelixClosingWeights(make_tuple(1, 7, 1, 1), rna2, &indices);
  expected = {mp(Features::HELIX_CLOSING_UA, 1)};
  EXPECT_EQ(expected, indices);

  HelixClosingWeights(make_tuple(4, 11, 1, 1), rna2, &indices);
  expected = {mp(Features::HELIX_CLOSING_UA, 1),
              mp(Features::HELIX_CLOSING_UC, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  HelixClosingWeights(make_tuple(4, 11, 1, 1), rna2, &indices);
  expected = {mp(Features::HELIX_CLOSING_UC, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  HelixClosingWeights(make_tuple(2, 6, 1, 1), rna2, &indices);
  expected = {mp(Features::HELIX_CLOSING_GG, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  HelixClosingWeights(make_tuple(1, 4, 1, 1), rna2, &indices);
  expected = {mp(Features::HELIX_CLOSING_UU, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  HelixClosingWeights(make_tuple(1, 5, 1, 1), rna2, &indices);
  expected = {mp(Features::HELIX_CLOSING_UG, 1)};
  EXPECT_EQ(expected, indices);
}

TEST_F(ModelTest, HelixExtendWeights) {
  Structure rna("............", "ACCACCTCCCGA", {});
  IndexVec indices;
  HelixExtendWeights(make_tuple(1, 1, 16, 1), rna, &indices);
  IndexVec expected = {mp(Features::EXTENSION_SCORE, 1)};
  EXPECT_EQ(expected, indices);
}

TEST_F(ModelTest, HelixStackingWeights) {
  Structure rna("..............", "AAGCAACGCCCCGA", {});
  IndexVec indices;

  HelixStackingWeights(make_tuple(1, 6, 2, 5), rna, &indices);
  IndexVec expected = {mp(Features::STACKING_PAIR_AAAA, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  HelixStackingWeights(make_tuple(3, 9, 4, 8), rna, &indices);
  expected = {mp(Features::STACKING_PAIR_AAAA + 109, 1)};
  EXPECT_EQ(expected, indices);
}

TEST_F(ModelTest, MultiMismatchWeights) {
  Structure rna("..............", "AACCAACCTCCCGA", {});
  IndexVec indices;

  MultiMismatchWeights(make_tuple(1, 6, 2, 5), rna, &indices);
  IndexVec expected = {mp(Features::SINGLE_BASE_PAIR_STACKING_LEFT_AAA, 1),
                       mp(Features::SINGLE_BASE_PAIR_STACKING_RIGHT_AAA, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  MultiMismatchWeights(make_tuple(1, 7, 2, 5), rna, &indices);
  expected = {mp(Features::SINGLE_BASE_PAIR_STACKING_LEFT_AAA + 4, 1),
              mp(Features::SINGLE_BASE_PAIR_STACKING_RIGHT_AAA + 4, 1)};
  EXPECT_EQ(expected, indices);

  // EXPECT_EQ(1, CalculateLeftOffset(make_tuple(0, 0, 1, IDX_NOT_SET)));

  // EXPECT_EQ(20, CalculateLeftOffset(make_tuple(1, 1, 0, IDX_NOT_SET)));
  // EXPECT_EQ(21, CalculateLeftOffset(make_tuple(1, 1, 1, IDX_NOT_SET)));
  // EXPECT_EQ(22, CalculateLeftOffset(make_tuple(1, 1, 2, IDX_NOT_SET)));

  // EXPECT_EQ(41, CalculateLeftOffset(make_tuple(2, 2, 1, IDX_NOT_SET)));
  // EXPECT_EQ(42, CalculateLeftOffset(make_tuple(2, 2, 2, IDX_NOT_SET)));
  // EXPECT_EQ(43, CalculateLeftOffset(make_tuple(2, 2, 3, IDX_NOT_SET)));

  // EXPECT_EQ(62, CalculateLeftOffset(make_tuple(3, 3, 2, IDX_NOT_SET)));
  // EXPECT_EQ(63, CalculateLeftOffset(make_tuple(3, 3, 3, IDX_NOT_SET)));

  // EXPECT_EQ(0, CalculateRightOffset(make_tuple(0, 0, IDX_NOT_SET, 0)));
  // EXPECT_EQ(1, CalculateRightOffset(make_tuple(0, 0, IDX_NOT_SET, 1)));

  // EXPECT_EQ(20, CalculateRightOffset(make_tuple(1, 1, IDX_NOT_SET, 0)));
  // EXPECT_EQ(21, CalculateRightOffset(make_tuple(1, 1, IDX_NOT_SET, 1)));
  // EXPECT_EQ(22, CalculateRightOffset(make_tuple(1, 1, IDX_NOT_SET, 2)));

  // EXPECT_EQ(41, CalculateRightOffset(make_tuple(2, 2, IDX_NOT_SET, 1)));
  // EXPECT_EQ(42, CalculateRightOffset(make_tuple(2, 2, IDX_NOT_SET, 2)));
  // EXPECT_EQ(43, CalculateRightOffset(make_tuple(2, 2, IDX_NOT_SET, 3)));

  // EXPECT_EQ(62, CalculateRightOffset(make_tuple(3, 3, IDX_NOT_SET, 2)));
  // EXPECT_EQ(63, CalculateRightOffset(make_tuple(3, 3, IDX_NOT_SET, 3)));
}

TEST_F(ModelTest, SeveralMultiWeights) {
  Structure rna("............", "ACCACCTCCCGA", {});
  IndexVec indices;

  MultiBaseWeights(make_tuple(1, 1, 1, 1), rna, &indices);
  IndexVec expected = {mp(Features::MULTI_BASE, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  MultiPairedWeights(make_tuple(1, 1, 16, 1), rna, &indices);
  expected = {mp(Features::MULTI_PAIRED, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  MultiUnpairedWeights(make_tuple(1, 1, 16, 1), rna, &indices);
  expected = {mp(Features::MULTI_UNPAIRED, 1)};
  EXPECT_EQ(expected, indices);

  indices.clear();
  MultiBaseWeights(make_tuple(1, 1, 16, 1), rna, &indices);
  MultiPairedWeights(make_tuple(1, 1, 16, 1), rna, &indices);
  MultiUnpairedWeights(make_tuple(1, 1, 16, 1), rna, &indices);
  expected = {mp(Features::MULTI_BASE, 1), mp(Features::MULTI_PAIRED, 1),
              mp(Features::MULTI_UNPAIRED, 1)};
  EXPECT_EQ(expected, indices);
}

// TEST_F(ModelTest, OuterUnpairedAndBranchWeights) {}
//
bool CompareIdxPair(pair<idx_t, score_t> a, pair<idx_t, score_t> b) {
  return a.first < b.first;
}

TEST_F(ModelTest, SingleWeights) {
  Structure rna("...................", "ACCATCTCCCGAACCAGCT", {});
  IndexVec indices;
  SingleWeights(make_tuple(3, 16, 5, 12), rna, &indices);
  //  i | j+1 | i+1 |  j | ip+1 | jp | ip | jp+1
  //  3 |  17 |   4 | 16 |    6 | 12 |  5 |   13
  IndexVec expected = {mp(Features::INTERNAL_LENGTH_1 + 5, 1),
                       mp(Features::INTERNAL_FULL_2_4, 1),
                       mp(Features::INTERNAL_ASYMMETRY_2, 1),
                       mp(Features::TERMINAL_MISMATCH_CGAA, 1),
                       mp(Features::TERMINAL_MISMATCH_CAUA, 1)};
  std::sort(expected.begin(), expected.end(), CompareIdxPair);
  std::sort(indices.begin(), indices.end(), CompareIdxPair);
  EXPECT_EQ(expected, indices);
}

TEST_F(ModelTest, InitializeConvertFeaturesNotAffected) {
  FeatureVec model;
  FeatureVec model2;
  for (int i = 0; i < 12; i++) {
    auto p = make_pair(i, drand48());
    model.insert(p);
    model2.insert(p);
  }
  InitializeFeatures(&model2);
  ConvertFeatures(&model2);
  for (auto it = model.cbegin(); it != model.cend(); ++it) {
    const idx_t feature_id = it->first;
    const score_t feature_value = it->second;
    // Add feature contribution to the score of this example.
    FeatureVec::const_iterator larger_entry = model2.find(feature_id);
    EXPECT_TRUE(larger_entry != model2.end());
    EXPECT_NEAR(larger_entry->second, feature_value, kEpsilon);
  }
}

// TEST_F(ModelTest, InitializeConvertFeatures) {
//  FeatureVec model;
//  FeatureVec model2;
//  for (int i = 0; i < 906; i++) {
//    auto p = make_pair(i, drand48());
//    model.insert(p);
//    model2.insert(p);
//  }
//  InitializeFeatures(&model2);
//  ConvertFeatures(&model2);
//  for (auto it = model.cbegin(); it != model.cend(); ++it) {
//    const idx_t feature_id = it->first;
//    const score_t feature_value = it->second;
//    // Add feature contribution to the score of this example.
//    FeatureVec::const_iterator larger_entry = model2.find(feature_id);
//    EXPECT_TRUE(larger_entry != model2.end());
//    EXPECT_NEAR(larger_entry->second, feature_value, kEpsilon);
//  }
//}

TEST_F(ModelTest, CalculateFeatures) {
  Structure rna("......", "ACCATC", {});
  FeatureVec features;
  CalculateFeatures(rna, &features);
  FeatureVec expected = {mp(0, 6)};
  EXPECT_EQ(expected, features);

  Structure rna2("(((..((...))....)))", "ACCATCTCCCGAACCAGCT", {});
  EXPECT_EQ(rna2.Size(), 20);
  CalculateFeatures(rna2, &features);
  expected = {};
  // EXPECT_EQ(expected, features);
}

}  // namespace
}  // namespace qrisp

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
