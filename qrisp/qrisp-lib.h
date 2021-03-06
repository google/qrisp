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

#ifndef QRISP_QRISP_LIB_H_
#define QRISP_QRISP_LIB_H_

#include "dataset-utils.h"
#include "qrisp/proto/config.pb.h"
#include "structure.h"
#include "utils.h"

using namespace std;

namespace qrisp {

void StartPrediction(const Dataset& dataset, const QRSPModel& model,
                     const ConfigMessage& config);

void StartTraining(const Dataset& training_set, const Dataset& holdout_set,
                   const QRSPModel& model, const ConfigMessage& config);

}  // namespace qrisp

#endif  // QRISP_QRISP_LIB_H_
