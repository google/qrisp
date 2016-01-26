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

#include "structure.h"
#include "qrisp/proto/structure.pb.h"
#include "utils.h"

#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <tuple>

namespace qrisp {
namespace {

// using testing::EqualsProto;
// using testing::InSequence;
// using testing::Pair;
// using testing::Return;
// using testing::UnorderedElementsAre;

// extern bool CheckForPseudoKnot(const Structure& structure);

class RNAToolsTest : public testing::Test {
 protected:
  // static string GetFilePath(const string& filename) {
  //  const string file_dir = file::JoinPath(
  //      FLAGS_test_srcdir, "experimental/users/fdb/qrsp/testdata");
  //    return (file::JoinPath(file_dir, filename));
  //}

  RNAToolsTest() {}
  ~RNAToolsTest() {}

  virtual void SetUp() {}
  virtual void TearDown() {}
};

static constexpr char kStructureAsciipb[] = R"(
  rows { pos: 1 }
  rows { pos: 2 }
  rows { pos: 3 }
  rows { pos: 4 }
  rows { pos: 5 }
  rows { pos: 6 }
  )";

TEST_F(RNAToolsTest, ParseProtoFromString) {
  Structure s;
  s.InitializeFromAsciiPb(kStructureAsciipb);
  // Structure t;
  // StructureMessage m;
  // m.ParseFromString(kStructureAsciipb);
  // LOG(INFO) << m.rows_size();
  // EXPECT_TRUE(t.InitializeFromProto(m));
}

TEST_F(RNAToolsTest, GetBasesAt) {
  {
    Structure rna(".(..).", "AGAACG", {});
    EXPECT_EQ(true, rna.IsValidStructure());
    EXPECT_EQ(std::make_tuple(0, 2, 0, 0),
              rna.GetBasesAt(std::make_tuple(1, 2, 3, 4)));
    EXPECT_EQ(std::make_tuple(0, 0, 1, 2),
              rna.GetBasesAt(std::make_tuple(3, 4, 5, 6)));
  }
  {
    Structure rna(".(..).", "AGAACG", {});
    EXPECT_EQ(true, rna.IsValidStructure());
    EXPECT_EQ(std::make_tuple(0, 2, 0, 0),
              rna.GetBasesAt(std::make_tuple(1, 2, 3, 4)));
    EXPECT_EQ(std::make_tuple(0, 0, 1, 2),
              rna.GetBasesAt(std::make_tuple(3, 4, 5, 6)));
  }
  {
    Structure rna(".(..).", "AGAACG", {});
    EXPECT_EQ(true, rna.IsValidStructure());
    EXPECT_EQ(std::make_tuple(0, 2, 0, 0),
              rna.GetBasesAt(std::make_tuple(1, 2, 3, 4)));
    EXPECT_EQ(std::make_tuple(0, 0, 1, 2),
              rna.GetBasesAt(std::make_tuple(3, 4, 5, 6)));
  }
  {
    Structure rna(".(..).", "AGAACG", {});
    EXPECT_EQ(true, rna.IsValidStructure());
    EXPECT_EQ(std::make_tuple(0, 2, 0, 0),
              rna.GetBasesAt(std::make_tuple(1, 2, 3, 4)));
    EXPECT_EQ(std::make_tuple(0, 0, 1, 2),
              rna.GetBasesAt(std::make_tuple(3, 4, 5, 6)));
  }
  {
    Structure rna(".(..).", "AGAACG", {});
    EXPECT_EQ(true, rna.IsValidStructure());
    EXPECT_EQ(std::make_tuple(0, 2, 0, 0),
              rna.GetBasesAt(std::make_tuple(1, 2, 3, 4)));
    EXPECT_EQ(std::make_tuple(0, 0, 1, 2),
              rna.GetBasesAt(std::make_tuple(3, 4, 5, 6)));
  }
}

// TEST_F(RNAToolsTest, ContainsPseudoKnot) {
//  {
//    Structure s("......", "AAAAAA", {});
//    EXPECT_EQ(false, s.ContainsPseudoKnot());
//  }
//  {
//    Structure s("((..(....)))", "AAAAAAAAAAAA", {});
//    EXPECT_EQ(false, s.ContainsPseudoKnot());
//  }
//  {
//    Structure s("((..)...(...)...).", "AAAAAAAAAAAAAAAAAA", {});
//    EXPECT_EQ(false, s.ContainsPseudoKnot());
//  }
//}

TEST_F(RNAToolsTest, IsValidStructure) {
  {
    Structure rna({}, "AGA", {});
    EXPECT_EQ(false, rna.IsValidStructure());
  }
  {
    Structure rna("(.)", "AGG", {});
    EXPECT_EQ(false, rna.IsValidStructure());
  }
  {
    Structure rna("(.....", "AGAACG", {});
    EXPECT_EQ(false, rna.IsValidStructure());
  }
  {
    Structure rna(".......", "ACTCACT", {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    EXPECT_EQ(true, rna.IsValidStructure());
  }
}

TEST_F(RNAToolsTest, CalculateSubstructure) {
  {
    Structure rna("(((..((...))....)))", "ACCATCTCCCGAACCAGCT", {});
    vector<Substructure> substructures;
    vector<bool> accessible_positions;
    rna.CalculateSubstructure(&substructures, &accessible_positions);
    vector<bool> expected_accessible_positions(accessible_positions.size(),
                                               true);
    expected_accessible_positions[0] = false;
    EXPECT_EQ(expected_accessible_positions, accessible_positions);
    EXPECT_EQ(5, substructures.size());
    pair<idx_t, idx_t> outer_pair = make_pair(1, 19);
    vector<pair<idx_t, idx_t>> inner_pairs = {make_pair(2, 18)};
    Substructure expected0(outer_pair, inner_pairs);
    EXPECT_EQ(expected0, substructures[0]);

    outer_pair = make_pair(2, 18);
    inner_pairs = {make_pair(3, 17)};
    Substructure expected1(outer_pair, inner_pairs);
    EXPECT_EQ(expected1, substructures[1]);

    outer_pair = make_pair(3, 17);
    inner_pairs = {make_pair(6, 12)};
    Substructure expected2(outer_pair, inner_pairs);
    EXPECT_EQ(expected2, substructures[2]);

    outer_pair = make_pair(6, 12);
    inner_pairs = {make_pair(7, 11)};
    Substructure expected3(outer_pair, inner_pairs);
    EXPECT_EQ(expected3, substructures[3]);

    outer_pair = make_pair(7, 11);
    inner_pairs = {};
    Substructure expected4(outer_pair, inner_pairs);
    EXPECT_EQ(expected4, substructures[4]);
  }
  {
    Structure rna("((((...))((....))))", "ACCATCTCCCGAACCAGCT", {});
    vector<Substructure> substructures;
    vector<bool> accessible_positions;
    rna.CalculateSubstructure(&substructures, &accessible_positions);
    vector<bool> expected_accessible_positions(accessible_positions.size(),
                                               true);
    expected_accessible_positions[0] = false;
    EXPECT_EQ(expected_accessible_positions, accessible_positions);
    EXPECT_EQ(6, substructures.size());
    pair<idx_t, idx_t> outer_pair = make_pair(1, 19);
    vector<pair<idx_t, idx_t>> inner_pairs = {make_pair(2, 18)};
    Substructure expected0(outer_pair, inner_pairs);
    EXPECT_EQ(expected0, substructures[0]);

    outer_pair = make_pair(2, 18);
    inner_pairs = {make_pair(3, 9), make_pair(10, 17)};
    Substructure expected1(outer_pair, inner_pairs);
    EXPECT_EQ(expected1, substructures[1]);

    outer_pair = make_pair(3, 9);
    inner_pairs = {make_pair(4, 8)};
    Substructure expected2(outer_pair, inner_pairs);
    EXPECT_EQ(expected2, substructures[2]);

    outer_pair = make_pair(4, 8);
    inner_pairs = {};
    Substructure expected3(outer_pair, inner_pairs);
    EXPECT_EQ(expected3, substructures[3]);
  }
}

// Now tests for auxiliary functions.
TEST_F(RNAToolsTest, PairingsToBrackets) {
  string brackets;
  PairingsToBrackets({IDX_NOT_SET, 0, 0, 0, 0, 0, 0}, &brackets);
  EXPECT_EQ("......", brackets);
  PairingsToBrackets({IDX_NOT_SET, 6, 0, 0, 0, 0, 1}, &brackets);
  EXPECT_EQ("(....)", brackets);
  PairingsToBrackets({IDX_NOT_SET, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                     &brackets);
  EXPECT_EQ("..............", brackets);
  PairingsToBrackets({IDX_NOT_SET, 14, 13, 0, 0, 12, 0, 0, 0, 0, 0, 0, 5, 2, 1},
                     &brackets);
  EXPECT_EQ("((..(......)))", brackets);
  PairingsToBrackets({IDX_NOT_SET, 7, 6, 0, 0, 0, 2, 1, 0, 0, 14, 0, 0, 0, 10},
                     &brackets);
  EXPECT_EQ("((...))..(...)", brackets);
}

TEST_F(RNAToolsTest, BracketsToPairings) {
  string brackets = "......";
  vector<idx_t> pairings;
  pairings.resize(brackets.size() + 1, IDX_NOT_SET);
  BracketsToPairings(brackets, &pairings);
  vector<idx_t> expected = {IDX_NOT_SET, 0, 0, 0, 0, 0, 0};
  EXPECT_EQ(expected, pairings);
  brackets = ".(..).";
  pairings.resize(brackets.size() + 1, IDX_NOT_SET);
  BracketsToPairings(brackets, &pairings);
  expected = {IDX_NOT_SET, 0, 5, 0, 0, 2, 0};
  EXPECT_EQ(expected, pairings);
}

TEST_F(RNAToolsTest, PairwiseLoss) {
  EXPECT_EQ(0, PairwiseLoss(Structure("......", "AAAAAA", {}),
                            Structure("......", "AAGGAA", {})));
  EXPECT_EQ(2, PairwiseLoss(Structure(".(..).", "AAAAAA", {}),
                            Structure("......", "AAGGAA", {})));
  EXPECT_EQ(4, PairwiseLoss(Structure("((..))", "AAAAAA", {}),
                            Structure("......", "AAGGAA", {})));
  EXPECT_EQ(3, PairwiseLoss(Structure(".(..).", "AAAAAA", {}),
                            Structure(".(...)", "AAGGAA", {})));
}

TEST_F(RNAToolsTest, HammingLoss) {
  EXPECT_EQ(0, HammingLoss(Structure("......", "AAAAAA", {}),
                           Structure("......", "AAGGAA", {})));
  EXPECT_EQ(2, HammingLoss(Structure(".(..).", "AAAAAA", {}),
                           Structure("......", "AAGGAA", {})));
  EXPECT_EQ(4, HammingLoss(Structure("((..))", "AAAAAA", {}),
                           Structure("......", "AAGGAA", {})));
  EXPECT_EQ(2, HammingLoss(Structure(".(..).", "AAAAAA", {}),
                           Structure(".(...)", "AAGGAA", {})));
}

}  // namespace
}  // namespace qrisp

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
