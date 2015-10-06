#include "dataset-utils.h"

#include <string>

namespace qrisp {

using namespace std;

// void LoadModelFromFile(const string& source, QRSPModel* model) {

// void SaveModelToFile(const string& destination, const QRSPModel& model) {

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

bool CompareKey(const pair<string, Structure>& s,
                const pair<string, Structure>& t) {
  return (s.first.compare(t.first) < 0);
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

bool IsComment(const string& line) {
  return (line.find("#") != string::npos || line.find("//") != string::npos);
}

void LoadSplit(const string& path, vector<string>* split) {
  // CHECK_OK(file::GetContents(source, &content, file::Defaults()));
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
  // map<string, RNASStructure> protomap;
  // SSTable_Utils::ReadUniqueProtoMapFromSSTableOrDie(source, &protomap);
  // for (const auto& elem : protomap) {
  //  data->insert(make_pair(elem.first, Structure(elem.second)));
  //}
}

void SaveDataset(const string& destination, const Dataset& data) {
  // map<string, RNASStructure> protomap;
  // for (const auto& elem : data) {
  //  RNASStructure structure;
  //  elem.second.ConvertToProto(&structure);
  //  protomap.insert(make_pair(elem.first, structure));
  //}
  // SSTable_Utils::WriteProtoMapToSSTableOrDie(protomap, destination,
  //                                           SSTableBuilderOptions());
}

std::tuple<double, double> SSMeasure(const Structure& ref,
                                     const Structure& pred) {
  const static auto null_tuple = std::make_tuple(0.0, 0.0);
  if (ref.GetSize() == 0) {
    return null_tuple;
  }
  CHECK(ref.GetSize() == pred.GetSize());

  int num_ref_pairs = 0;
  int num_pred_pairs = 0;
  int correctly_paired = 0;

  for (int i = 1; i < pred.GetSize(); i++) {
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
