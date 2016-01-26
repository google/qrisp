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

#ifndef QRISP_CLUSTER_H_
#define QRISP_CLUSTER_H_

#include "utils.h"

#include <map>
#include <utility>
#include <vector>

namespace qrisp {

using namespace std;

typedef map<pair<int, int>, int> PairHistogram;

// Each base pair comes with a frequency,
// (i, i, frequency).
bool CalculateCentroid(const vector<vector<idx_t>> pairings_multiset,
                       vector<idx_t>* centroid);

}  // namespace qrisp

#endif  // QRISP_CLUSTER_H_
