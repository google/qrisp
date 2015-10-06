#ifndef QRISP_PERFORMANCE_H_
#define QRISP_PERFORMANCE_H_

#include "proto/config.pb.h"
#include "utils.h"

namespace qrisp {

typedef vector<vector<vector<int>>> Predictions;
typedef map<int, pair<double, double>> PerformanceMetrics;

bool EstimatePerformanceOnHoldout(const Dataset& data,
                                  const Predictions& predictions,
                                  const vector<int> cluster_sizes,
                                  PerformanceMetrics* performance);
}  // namespace qrisp
#endif  // EXPERIMENTAL_USERS_FDB_QRSP_PERFORMANCE_H_
