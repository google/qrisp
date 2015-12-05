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
#include "utils.h"

#include <gtest/gtest.h>

namespace qrisp {
namespace {

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
  EXPECT_NEAR(0.0, limits[0], epsilon);
  EXPECT_NEAR(3.333333, limits[1], epsilon);
  EXPECT_NEAR(6.666666, limits[2], epsilon);
  EXPECT_NEAR(10.0, limits[3], epsilon);

  linspace(-1.0, 1.0, 3, &limits);
  EXPECT_NEAR(-1.0, limits[0], epsilon);
  EXPECT_NEAR(0.0, limits[1], epsilon);
  EXPECT_NEAR(1.0, limits[2], epsilon);

  // TODO(fabio): Add test for logspace as well.
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