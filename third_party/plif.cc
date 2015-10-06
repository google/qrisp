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

#include <tuple>

namespace qrisp {

void EvaluatePlifAt(const ScoreVec& limits, const double& x, const idx_t offset,
                    IndexVec* indices) {
  const size_t len = limits.size();
  if (x < limits[0]) {
    indices->push_back(make_pair(0, 1.0));
    return;
  }
  if (x >= limits[len - 1]) {
    indices->push_back(make_pair(len - 1, 1.0));
    return;
  }
  for (size_t i = 0; i < len - 1; i++) {
    double limits_i = limits[i];
    double limits_ip = limits[i + 1];
    if (limits_i <= x && x < limits_ip) {
      double denom = limits_ip - limits_i;
      indices->push_back(make_pair(offset + i, (limits_ip - x) / denom));
      indices->push_back(make_pair(offset + i + 1, (x - limits_i) / denom));
    }
  }
}

typedef std::tuple<idx_t, double, double, double> coord_t;

coord_t EvaluateOneDimension(const ScoreVec& limits, const double& x) {
  const size_t len = limits.size();
  if (x < limits[0]) {
    return std::make_tuple(0, 1.0, 0.0, 1.0);
  }
  if (x >= limits[len - 1]) {
    // TODO(fdb): was -2 here but should be -1 I guess.
    return std::make_tuple(len - 2, 1.0, 1.0, 0.0);
  }
  for (size_t i = 0; i < len - 1; i++) {
    double limits_i = limits[i];
    double limits_ip = limits[i + 1];
    if (limits_i <= x && x < limits_ip) {
      return std::make_tuple(i, limits_ip - limits_i, x - limits_i,
                             limits_ip - x);
    }
  }
  CHECK(false);
  coord_t inter(0, 0.0, 0.0, 0.0);
  return inter;
}

void Evaluate2DPlifAt(const ScoreVec& limits1, const ScoreVec& limits2,
                      const double& x, const double& y, const idx_t offset,
                      IndexVec* indices) {
  const size_t len1 = limits1.size();
  double dx, dx1, dx2;
  double dy, dy1, dy2;
  idx_t i, j;
  std::tie(i, dx, dx1, dx2) = EvaluateOneDimension(limits1, x);
  std::tie(j, dy, dy1, dy2) = EvaluateOneDimension(limits2, y);
  const double denom = dx * dy;
  indices->push_back(make_pair(offset + j * len1 + i, (dx2 * dy2) / denom));
  indices->push_back(make_pair(offset + j * len1 + i + 1, (dx1 * dy2) / denom));
  indices->push_back(
      make_pair(offset + (j + 1) * len1 + i + 1, (dx1 * dy1) / denom));
  indices->push_back(
      make_pair(offset + (j + 1) * len1 + i, (dx2 * dy1) / denom));
}

}  // namespace qrisp
