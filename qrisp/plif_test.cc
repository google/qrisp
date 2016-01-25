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

#include "plif.h"
#include <gtest/gtest.h>
#include "utils.h"

namespace qrisp {
namespace {

constexpr const float kEpsilon = 1e-6;

class PlifTest : public testing::Test {
 protected:
  PlifTest() {}
  ~PlifTest() {}

  virtual void SetUp() {}
  virtual void TearDown() {}
};

TEST_F(PlifTest, AuxiliaryFunctions) {
  ScoreVec limits;
  linspace(0.0, 10.0, 4, &limits);
  EXPECT_NEAR(0.0, limits[0], kEpsilon);
  EXPECT_NEAR(3.333333, limits[1], kEpsilon);
  EXPECT_NEAR(6.666666, limits[2], kEpsilon);
  EXPECT_NEAR(10.0, limits[3], kEpsilon);

  linspace(-1.0, 1.0, 3, &limits);
  EXPECT_NEAR(-1.0, limits[0], kEpsilon);
  EXPECT_NEAR(0.0, limits[1], kEpsilon);
  EXPECT_NEAR(1.0, limits[2], kEpsilon);
}

TEST_F(PlifTest, EvaluatePlifAt) {
  const ScoreVec limits = {0.0, 1.0, 2.0, 3.0, 4.0};
  IndexVec indices;
  EvaluatePlifAt(limits, 0.5, 0, &indices);
  const IndexVec expected = {mp(0, 0.5), mp(1, 0.5)};
  EXPECT_EQ(expected, indices);
  const IndexVec expected2 = {mp(0, 0.9), mp(1, 0.1)};
  indices.clear();
  EvaluatePlifAt(limits, 0.1, 0, &indices);
  EXPECT_EQ(expected2, indices);
}

}  // namespace
}  // namespace qrisp

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
