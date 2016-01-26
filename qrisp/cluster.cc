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

#include "cluster.h"

namespace qrisp {
bool CalculateCentroid(const vector<vector<idx_t>> pairings_multiset,
                       vector<idx_t>* centroid) {
  if (pairings_multiset.size() == 0) {
    return false;
  }
  const int fixed_size = pairings_multiset[0].size();
  centroid->clear();
  centroid->assign(fixed_size, 0);
  (*centroid)[0] = -1;
  if (fixed_size != centroid->size()) {
    return false;
  }
  if (pairings_multiset.size() == 1) {
    *centroid = pairings_multiset[0];
    return true;
  }
  PairHistogram pair_frequencies;
  for (const auto& pairings : pairings_multiset) {
    if (pairings.size() != fixed_size) {
      return false;
    }
    for (size_t i = 1; i < pairings.size(); i++) {
      if (pairings[i] != 0 && i < pairings[i]) {
        auto current_pair = make_pair(i, pairings[i]);
        if (pair_frequencies.find(current_pair) != pair_frequencies.end()) {
          pair_frequencies[current_pair]++;
        } else {
          pair_frequencies[current_pair] = 1;
        }
      }
    }
  }
  const int num_structures = pairings_multiset.size();
  const int threshold = num_structures / 2;
  for (const auto& p : pair_frequencies) {
    if (p.second > threshold) {
      (*centroid)[p.first.first] = p.first.second;
      (*centroid)[p.first.second] = p.first.first;
    }
  }
  return true;
}

}  // namespace qrisp
