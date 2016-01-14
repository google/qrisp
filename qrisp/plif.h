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

/*
   3 * 30 supporting points: x_1 ... x_30 => 3 * 30 parameters (1 .. 90):
     y_1 ... y_30
   piecewise linear function:
   Take score from SVM vector (don_supp: case 1, acc_supp: case 2) and compute
   length of intron: case 3 these are our values x

          | y_1                                          if x <= x_1
          |
          |        x_i+1 - x              x - x_i
   f(x) = | y_i * -----------  +  y_i+1 * -----------    if x_i <= x <= x_i+1
          |       x_i+1 - x_i             x_i+1 - x_i
          |
          | y_30                                         if x_n <= x

    y_i and y_i+1 parameters, so the fractions are saved in the weight vectors!
 *
 */

#ifndef QRISP_PLIF_H__
#define QRISP_PLIF_H__

#include "utils.h"

#include <cmath>
#include <vector>

namespace qrisp {

// Auxiliary function that helps creating a vector of support points for
// piecewise-linear functions. Given a lower and an upper bound together with
// the number of buckets this function creates a vector that contains values
// that separate the interval [lower_bound, upper_bound] into num_buckets
// buckets of equal size.
inline void linspace(const double& lower_bound, const double& upper_bound,
                     const int& num_buckets, ScoreVec* limits) {
  limits->resize(num_buckets, 0.0);
  (*limits)[0] = lower_bound;
  (*limits)[num_buckets - 1] = upper_bound;
  const double step_size = (upper_bound - lower_bound) / (num_buckets - 1);
  for (int i = 1; i < num_buckets - 1; i++)
    (*limits)[i] = lower_bound + (i * step_size);
}

// This function can be used to calculate the update values for a 1d plif
// on-the-fly without using the Plif class above.
void EvaluatePlifAt(const ScoreVec& limits, const double& x, const idx_t offset,
                    IndexVec* indices);

// This function can be used to calculate the update values for a 2d plif
// on-the-fly without using the Plif2D class above.
//
//
//     +----------------------------+
//     |                  |         |
//     |                  |         |
//     |                  |         |
//     |------------------+---------|
//     |                  |         |
//     |                  |         |
//     |                  |         |
//     |                  |         |
//     |                  |         |
//     |                  |         |
//   Y |                  |         |
//     |                  |         |
//     +----------------------------+
//        X
//
//
//
void Evaluate2DPlifAt(const ScoreVec& limits1, const ScoreVec& limits2,
                      const double& x, const double& y, const idx_t offset,
                      IndexVec* indices);

}  // namespace qrisp

#endif  // QRISP_PLIF_H__
