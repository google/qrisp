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

#include "utils.h"

namespace qrisp {

// bool IsValidConfig(const bool& training, const QRFPConfig& config) {
//  if(training) {
//    // Check for settings needed during training phase.
//    if (config.result_dir().empty()) {
//      cout << "Error: You have to supply a result directory!" << endl;
//      exit(EXIT_FAILURE);
//    }
//    if(config.training_algorithm() == "mosek-qp" ||
//       config.training_algorithm() == "subgradient") {
//      if(!(config.rescaling_method() == "margin" ||
//           config.rescaling_method() == "slack" ||
//           config.rescaling_method() == "none")) {
//        cout << "Error: When using subgradient or mosek-qp algorithms "
//             "rescaling_method has to be one of:" << endl;
//        cout << "margin|slack|none" << endl;
//        exit(EXIT_FAILURE);
//      }
//    }
//    if (config.score_lower_bound() < -10000 || config.score_lower_bound() >
//    10000) {
//      cout << "Error: The lower bound for the PARS score should range between
//      "
//              "-10000 and +10000." << endl;
//      exit(EXIT_FAILURE);
//    }
//    if (config.score_upper_bound() < -10000 || config.score_upper_bound() >
//    10000) {
//      cout << "Error: The upper bound for the PARS score should range between
//      "
//              "-10000 and +10000." << endl;
//      exit(EXIT_FAILURE);
//    }
//  } else {
//    // Check for settings needed during prediction phase.
//  }
//  return true;
//}
//
void SplitStringToInts(const string& input, vector<int>* output) {
  // TODO(fdb): Implement.
  // Fetch cluster sizes.
  // char* p = strtok(input.c_str(), ",");
  // double value;
  // while (p != NULL) {
  //  p = strtok(NULL, ","
  //  std::stod(string(p), &value);
  //  output->push_back(static_cast<int>(value));
  //}
}

}  // namespace qrisp
