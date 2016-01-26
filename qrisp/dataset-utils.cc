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

#include "dataset-utils.h"

namespace qrisp {

namespace {

bool IsComment(const string& line) {
  return (line.find("#") != string::npos || line.find("//") != string::npos);
}

}  // namespace

void InitializeFeatureVecFromModel(const QRSPModel& model, FeatureVec* fv) {
  for (const auto& feature : model.feature()) {
    auto insert_ret = fv->insert(make_pair(feature.id(), feature.value()));
    if (!insert_ret.second) {
      LOG(ERROR) << "Feature id appears twice!";
    }
  }
}

void InitializeModelFromFeatureVec(const FeatureVec& fv, QRSPModel* model) {
  for (const auto& feature : fv) {
    QRSPModel::Feature* f = model->mutable_feature()->Add();
    f->set_id(feature.first);
    f->set_value(feature.second);
  }
}

void SubselectData(const vector<string>& split, const Dataset& input,
                   Dataset* output) {
  for (const auto& id : split) {
    const auto& elem = input.find(id);
    if (elem != input.end()) {
      output->insert(*elem);
    }
  }
}

void LoadSplit(const string& path, vector<string>* split) {
  std::fstream fs;
  fs.open(path.c_str(), std::fstream::in);
  if (!fs.is_open()) {
    LOG(ERROR) << "Could not open file: " << path;
  }
  string line;
  while (getline(fs, line)) {
    if (IsComment(line)) {
      continue;
    }
    split->push_back(line);
  }
  fs.close();
}

void LoadDataset(const string& source, Dataset* data) {
  StructureSetMessage structures;
  LoadProtoFromFile(source, &structures);
  for (const auto& elem : structures.item()) {
    auto it = data->insert(make_pair(elem.id(), Structure()));
    if (!it.second) {
      LOG(ERROR) << "Element not inserted: " << elem.id();
      continue;
    }
    it.first->second.InitializeFromProto(elem.structure());
  }
}

void SaveDataset(const string& destination, const Dataset& data) {
  StructureSetMessage structures;
  for (const auto& elem : data) {
    const auto item = structures.add_item();
    item->set_id(elem.first);
    elem.second.ConvertToProto(item->mutable_structure());
  }
}

std::tuple<double, double> SSMeasure(const Structure& ref,
                                     const Structure& pred) {
  const static auto null_tuple = std::make_tuple(0.0, 0.0);
  if (ref.Size() == 0) {
    return null_tuple;
  }
  CHECK(ref.Size() == pred.Size());

  int num_ref_pairs = 0;
  int num_pred_pairs = 0;
  int correctly_paired = 0;

  for (int i = 1; i < pred.Size(); i++) {
    if (pred(i) != 0) {
      num_pred_pairs++;
    }
    if (ref(i) != 0) {
      num_ref_pairs++;
    }
    if (pred(i) != 0 && pred(i) == ref(i)) {
      correctly_paired++;
    }
  }
  return std::make_tuple(
      (num_ref_pairs == 0) ? 1 : (double)correctly_paired / num_ref_pairs,
      (num_pred_pairs == 0) ? 1 : (double)correctly_paired / num_pred_pairs);
}

}  // namespace qrisp
