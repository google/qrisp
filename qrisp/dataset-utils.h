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

#ifndef QRISP_DATASET_UTILS_H_
#define QRISP_DATASET_UTILS_H_

#include "qrisp/proto/config.pb.h"
#include "qrisp/proto/parameters.pb.h"
#include "qrisp/proto/structure.pb.h"
#include "structure.h"
#include "utils.h"

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>
#include <algorithm>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

namespace qrisp {
// Load a protocol buffer message stored in ASCII format.
template <typename T>
void LoadProtoFromFile(const string& source, T* proto) {
  std::ifstream input_stream;
  input_stream.open(source.c_str(), std::ios::in);
  if (!input_stream.is_open()) {
    LOG(ERROR) << "Failed to open input stream for file:\n" << source;
  }
  google::protobuf::io::IstreamInputStream file_input(&input_stream);
  if (!google::protobuf::TextFormat::Parse(&file_input, proto)) {
    LOG(ERROR) << "Failed to parse proto buffer from file:\n" << source;
  }
  input_stream.close();
}

// Dump the content of a protocol buffer message to a text file (proto ASCII
// format).
template <typename T>
void SaveProtoToFile(const string& destination, const T& proto) {
  std::ofstream output_stream;
  output_stream.open(destination.c_str(), std::ios::out | std::ios::binary);
  if (!output_stream.is_open()) {
    LOG(ERROR) << "Failed to open output stream for file:\n" << destination;
  }
  google::protobuf::io::OstreamOutputStream file_output(&output_stream);
  if (!google::protobuf::TextFormat::Print(proto, &file_output)) {
    LOG(ERROR) << "Failed to serialize proto buffer to stream:\n"
               << proto.DebugString();
  }
}

void LoadSplit(const string& source, vector<string>* split);

void LoadDataset(const string& source, Dataset* data);

void InitializeFeatureVecFromModel(const QRSPModel& model, FeatureVec* fv);

void InitializeModelFromFeatureVec(const FeatureVec& fv, QRSPModel* model);

void SubselectData(const vector<string>& split, const Dataset& input,
                   Dataset* output);

void SaveDataset(const string& destination, const Dataset& data);

std::tuple<double, double> SSMeasure(const Structure& ref,
                                     const Structure& pred);

}  // namespace qrisp

#endif  // QRISP_DATASET_UTILS_H_
