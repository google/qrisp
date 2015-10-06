#ifndef QRISP_CLUSTER_H_
#define QRISP_CLUSTER_H_

#include <map>
#include <utility>
#include <vector>

namespace qrisp {

using namespace std;

typedef map<pair<int, int>, int> PairHistogram;

// Each base pair comes with a frequency,
// (i, i, frequency).
bool CalculateCentroid(const vector<vector<int>> pairings_multiset,
                       vector<int>* centroid);

}  // namespace qrisp

#endif  // QRISP_CLUSTER_H_
