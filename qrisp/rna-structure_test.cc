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

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <tuple>

#include <gtest/gtest.h>

#include "qrisp/proto/structure.pb.h"
#include "rna-structure.h"
#include "utils.h"

namespace qrisp {
namespace {

//using testing::EqualsProto;
//using testing::InSequence;
//using testing::Pair;
//using testing::Return;
//using testing::UnorderedElementsAre;

//extern bool CheckForPseudoKnot(const Structure& structure);

class RNAToolsTest : public testing::Test {
 protected:
  //static string GetFilePath(const string& filename) {
  //  const string file_dir = file::JoinPath(
  //      FLAGS_test_srcdir, "experimental/users/fdb/qrsp/testdata");
  //    return (file::JoinPath(file_dir, filename));
  //}

  RNAToolsTest() {}
  ~RNAToolsTest() {}

  virtual void SetUp() {}
  virtual void TearDown() {}
};

constexpr const char kStructureAsciipb[] = R"(
  length: 6
  base: 0
  base: 0
  base: 0
  base: 0
  base: 0
  base: 0
)";

TEST_F(RNAToolsTest, Initialization) {
  Structure s;
  EXPECT_TRUE(s.Initialize("......", "AAAAAA", {}));
  Structure t;
  StructureMessage m;
  m.ParseFromString(kStructureAsciipb);
  //EXPECT_TRUE(t.InitializeFromProto(m));

}

//TEST_F(RNAToolsTest, LoadTestdata) {
//  StructureMessages expected_results;
//  ReadDataFromFile<StructureMessages>("structures.txt", &expected_results);
//  EXPECT_EQ(expected_results.instance_size(), 1);
//}
//k
//kTEST_F(RNAToolsTest, Constructors) {
//k  {
//k  Structure rna_a("......", "AAAAAA", {});
//k  vector<idx_t> pairings = {IDX_NOT_SET, 0, 0, 0, 0, 0, 0};
//k  Structure rna_b(pairings, "AAAAAA", {});
//k  EXPECT_EQ(rna_a, rna_b);
//k  }
//k  {
//k  Structure rna_a(".(..).", "AGAACA", {});
//k  vector<idx_t> pairings = {IDX_NOT_SET, 0, 5, 0, 0, 2, 0};
//k  Structure rna_b(pairings, "AGAACA", {});
//k  }
//k  {
//k  Structure rna_a(".(..).", "AGAACA", { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 });
//k  }
//k}
//k
//kTEST_F(RNAToolsTest, GetBasesAt) {
//k  {
//k    Structure rna(".(..).", "AGAACG", {});
//k    EXPECT_EQ(true, rna.IsValidStructure());
//k    EXPECT_EQ(std::make_tuple(0, 2, 0, 0),
//k              rna.GetBasesAt(std::make_tuple(1, 2, 3, 4)));
//k    EXPECT_EQ(std::make_tuple(0, 0, 1, 2),
//k              rna.GetBasesAt(std::make_tuple(3, 4, 5, 6)));
//k  }
//k  {
//k    Structure rna(".(..).", "AGAACG", {});
//k    EXPECT_EQ(true, rna.IsValidStructure());
//k    EXPECT_EQ(std::make_tuple(0, 2, 0, 0),
//k              rna.GetBasesAt(std::make_tuple(1, 2, 3, 4)));
//k    EXPECT_EQ(std::make_tuple(0, 0, 1, 2),
//k              rna.GetBasesAt(std::make_tuple(3, 4, 5, 6)));
//k  }
//k  {
//k    Structure rna(".(..).", "AGAACG", {});
//k    EXPECT_EQ(true, rna.IsValidStructure());
//k    EXPECT_EQ(std::make_tuple(0, 2, 0, 0),
//k              rna.GetBasesAt(std::make_tuple(1, 2, 3, 4)));
//k    EXPECT_EQ(std::make_tuple(0, 0, 1, 2),
//k              rna.GetBasesAt(std::make_tuple(3, 4, 5, 6)));
//k  }
//k  {
//k    Structure rna(".(..).", "AGAACG", {});
//k    EXPECT_EQ(true, rna.IsValidStructure());
//k    EXPECT_EQ(std::make_tuple(0, 2, 0, 0),
//k              rna.GetBasesAt(std::make_tuple(1, 2, 3, 4)));
//k    EXPECT_EQ(std::make_tuple(0, 0, 1, 2),
//k              rna.GetBasesAt(std::make_tuple(3, 4, 5, 6)));
//k  }
//k  {
//k    Structure rna(".(..).", "AGAACG", {});
//k    EXPECT_EQ(true, rna.IsValidStructure());
//k    EXPECT_EQ(std::make_tuple(0, 2, 0, 0),
//k              rna.GetBasesAt(std::make_tuple(1, 2, 3, 4)));
//k    EXPECT_EQ(std::make_tuple(0, 0, 1, 2),
//k              rna.GetBasesAt(std::make_tuple(3, 4, 5, 6)));
//k  }
//k}
//k
//kTEST_F(RNAToolsTest, ContainsPseudoKnot) {
//k  {
//k  Structure s("......", "AAAAAA", {});
//k  EXPECT_EQ(false, s.ContainsPseudoKnot());
//k  }
//k  //{
//k  //Structure s("((..(....)))", "AAAAAAAAAAAA", {});
//k  //EXPECT_EQ(false, s.ContainsPseudoKnot());
//k  //}
//k  //{
//k  //Structure s("((..)...(...)...).", "AAAAAAAAAAAAAAAAAA", {});
//k  //EXPECT_EQ(false, s.ContainsPseudoKnot());
//k  //}
//k}
//k
//k//TEST_F(RNAToolsTest, IsValidStructure) {
//k//  {
//k//    Structure rna({}, "AGA", {});
//k//    EXPECT_EQ(false, rna.IsValidStructure());
//k//  }
//k//  {
//k//    Structure rna("(.)", "AGG", {});
//k//    EXPECT_EQ(false, rna.IsValidStructure());
//k//  }
//k//  {
//k//    Structure rna("(.....", "AGAACG", {});
//k//    EXPECT_EQ(false, rna.IsValidStructure());
//k//  }
//k//}
//k
//kTEST_F(RNAToolsTest, CalculateSubstructure) {
//k  {
//k  Structure rna("(((..((...))....)))", "ACCATCTCCCGAACCAGCT", {});
//k  vector<Substructure> substructures;
//k  vector<bool> accessible_positions;
//k  rna.CalculateSubstructure(&substructures, &accessible_positions);
//k  vector<bool> expected_accessible_positions(accessible_positions.size(), true);
//k  expected_accessible_positions[0] = false;
//k  EXPECT_EQ(expected_accessible_positions, accessible_positions);
//k  EXPECT_EQ(5, substructures.size());
//k  pair<idx_t, idx_t> outer_pair = make_pair(1, 19);
//k  vector<pair<idx_t, idx_t>> inner_pairs = {make_pair(2, 18)};
//k  Substructure expected0(outer_pair, inner_pairs);
//k  EXPECT_EQ(expected0, substructures[0]);
//k
//k  outer_pair = make_pair(2, 18);
//k  inner_pairs = {make_pair(3, 17)};
//k  Substructure expected1(outer_pair, inner_pairs);
//k  EXPECT_EQ(expected1, substructures[1]);
//k
//k  outer_pair = make_pair(3, 17);
//k  inner_pairs = {make_pair(6, 12)};
//k  Substructure expected2(outer_pair, inner_pairs);
//k  EXPECT_EQ(expected2, substructures[2]);
//k
//k  outer_pair = make_pair(6, 12);
//k  inner_pairs = {make_pair(7, 11)};
//k  Substructure expected3(outer_pair, inner_pairs);
//k  EXPECT_EQ(expected3, substructures[3]);
//k
//k  outer_pair = make_pair(7, 11);
//k  inner_pairs = {};
//k  Substructure expected4(outer_pair, inner_pairs);
//k  EXPECT_EQ(expected4, substructures[4]);
//k  }
//k  {
//k  Structure rna("((((...))((....))))", "ACCATCTCCCGAACCAGCT", {});
//k  vector<Substructure> substructures;
//k  vector<bool> accessible_positions;
//k  rna.CalculateSubstructure(&substructures, &accessible_positions);
//k  vector<bool> expected_accessible_positions(accessible_positions.size(), true);
//k  expected_accessible_positions[0] = false;
//k  EXPECT_EQ(expected_accessible_positions, accessible_positions);
//k  EXPECT_EQ(6, substructures.size());
//k  pair<idx_t, idx_t> outer_pair = make_pair(1, 19);
//k  vector<pair<idx_t, idx_t>> inner_pairs = {make_pair(2, 18)};
//k  Substructure expected0(outer_pair, inner_pairs);
//k  EXPECT_EQ(expected0, substructures[0]);
//k
//k  outer_pair = make_pair(2, 18);
//k  inner_pairs = {make_pair(3, 9), make_pair(10, 17)};
//k  Substructure expected1(outer_pair, inner_pairs);
//k  EXPECT_EQ(expected1, substructures[1]);
//k
//k  outer_pair = make_pair(3, 9);
//k  inner_pairs = {make_pair(4, 8)};
//k  Substructure expected2(outer_pair, inner_pairs);
//k  EXPECT_EQ(expected2, substructures[2]);
//k
//k  outer_pair = make_pair(4, 8);
//k  inner_pairs = {};
//k  Substructure expected3(outer_pair, inner_pairs);
//k  EXPECT_EQ(expected3, substructures[3]);
//k  }
//k}
//k
//k// Now tests for auxiliary functions.
//k
//kTEST_F(RNAToolsTest, PairingsToBrackets) {
//k  string brackets;
//k  PairingsToBrackets({ IDX_NOT_SET, 0, 0, 0, 0, 0, 0 }, &brackets);
//k  EXPECT_EQ("......", brackets);
//k  PairingsToBrackets({ IDX_NOT_SET, 6, 0, 0, 0, 0, 1 }, &brackets);
//k  EXPECT_EQ("(....)", brackets);
//k  PairingsToBrackets({ IDX_NOT_SET, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
//k                     &brackets);
//k  EXPECT_EQ("..............", brackets);
//k  PairingsToBrackets({ IDX_NOT_SET, 14, 13, 0, 0, 12, 0, 0, 0, 0, 0, 0, 5, 2,
//k                       1 },
//k                     &brackets);
//k  EXPECT_EQ("((..(......)))", brackets);
//k  PairingsToBrackets({ IDX_NOT_SET, 7, 6, 0, 0, 0, 2, 1, 0, 0, 14, 0, 0, 0,
//k                       10 },
//k                     &brackets);
//k  EXPECT_EQ("((...))..(...)", brackets);
//k}
//k
//kTEST_F(RNAToolsTest, BracketsToPairings) {
//k  string brackets = "......";
//k  vector<idx_t> pairings;
//k  pairings.resize(brackets.size() + 1, IDX_NOT_SET);
//k  BracketsToPairings(brackets, &pairings);
//k  vector<idx_t> expected = { IDX_NOT_SET, 0, 0, 0, 0, 0, 0 };
//k  EXPECT_EQ(expected, pairings);
//k  brackets = ".(..).";
//k  pairings.resize(brackets.size() + 1, IDX_NOT_SET);
//k  BracketsToPairings(brackets, &pairings);
//k  expected = { IDX_NOT_SET, 0, 5, 0, 0, 2, 0 };
//k  EXPECT_EQ(expected, pairings);
//k}

//TEST_F(RNAToolsTest, PairwiseLoss) {
//  EXPECT_EQ(0, PairwiseLoss(Structure("......", "AAAAAA", {}),
//                            Structure("......", "AAGGAA", {})));
//  EXPECT_EQ(2, PairwiseLoss(Structure(".(..).", "AAAAAA", {}),
//                            Structure("......", "AAGGAA", {})));
//  EXPECT_EQ(4, PairwiseLoss(Structure("((..))", "AAAAAA", {}),
//                            Structure("......", "AAGGAA", {})));
//  EXPECT_EQ(3, PairwiseLoss(Structure(".(..).", "AAAAAA", {}),
//                            Structure(".(...)", "AAGGAA", {})));
//}
//
//TEST_F(RNAToolsTest, HammingLoss) {
//  EXPECT_EQ(0, HammingLoss(Structure("......", "AAAAAA", {}),
//                           Structure("......", "AAGGAA", {})));
//  EXPECT_EQ(2, HammingLoss(Structure(".(..).", "AAAAAA", {}),
//                           Structure("......", "AAGGAA", {})));
//  EXPECT_EQ(4, HammingLoss(Structure("((..))", "AAAAAA", {}),
//                           Structure("......", "AAGGAA", {})));
//  EXPECT_EQ(2, HammingLoss(Structure(".(..).", "AAAAAA", {}),
//                           Structure(".(...)", "AAGGAA", {})));
//}

}  // namespace
}  // namespace qrisp

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
