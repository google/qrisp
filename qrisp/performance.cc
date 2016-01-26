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

#include <numeric>

#include "cluster.h"
#include "dataset-utils.h"
#include "performance.h"
#include "structure.h"

namespace qrisp {

bool EstimatePerformanceOnHoldout(const Dataset& data,
                                  const Predictions& predictions,
                                  const vector<idx_t> cluster_sizes,
                                  PerformanceMetrics* performance) {
  if (data.size() == 0 || cluster_sizes.size() == 0 ||
      data.size() != predictions.size()) {
    return false;
  }
  const int num_distinct_clusters = cluster_sizes.size();
  vector<vector<double>> all_sensitivity(num_distinct_clusters,
                                         vector<double>());
  vector<vector<double>> all_specificity(num_distinct_clusters,
                                         vector<double>());
  int j = 0;
  for (const auto& elem : data) {
    const Structure& input = elem.second;
    const auto& pairings = predictions[j];
    // Estimate values for different centroid cluster sizes.
    for (int i = 0; i < num_distinct_clusters; i++) {
      const int n = cluster_sizes[i];
      if (n <= pairings.size()) {
        vector<vector<idx_t>> best_n(pairings.begin(), pairings.begin() + n);
        // for (const auto& elem : best_n) {
        //  string out;
        //  PairingsToBrackets(elem, &out);
        //  LOG(INFO) << out;
        //}
        vector<idx_t> centroid;
        CHECK(CalculateCentroid(best_n, &centroid))
            << "Error: Could not calc. centroid for instance " << elem.first;
        string centroid_brac;
        PairingsToBrackets(centroid, &centroid_brac);
        Structure decoded_structure(input);
        decoded_structure.SetPairings(centroid);
        double sens = 0.0;
        double spec = 0.0;
        std::tie(sens, spec) = SSMeasure(input, decoded_structure);
        all_sensitivity[i].push_back(sens);
        all_specificity[i].push_back(spec);
      }
    }
    j++;
  }
  for (int i = 0; i < num_distinct_clusters; i++) {
    const auto& sensitivity = all_sensitivity[i];
    const auto& specificity = all_specificity[i];
    double sens = 0.0;
    double spec = 0.0;
    LOG(INFO) << "Performances: ";
    if (sensitivity.size() > 0 && specificity.size() > 0) {
      sens = std::accumulate(sensitivity.begin(), sensitivity.end(), 0.0) /
             static_cast<double>(sensitivity.size());
      spec = std::accumulate(specificity.begin(), specificity.end(), 0.0) /
             static_cast<double>(specificity.size());
      for (int i = 0; i < sensitivity.size(); i++) {
        LOG(INFO) << sensitivity[i] << " " << specificity[i];
      }
    }
    const int n = cluster_sizes[i];
    LOG(INFO) << "Holdout sensitivity / specificity @" << n << ": " << sens
              << " / " << spec;
    performance->emplace(std::make_pair(n, std::make_pair(sens, spec)));
  }
  return true;
}

}  // namespace qrisp
