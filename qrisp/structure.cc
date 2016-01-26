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

#include "structure.h"

#include <glog/logging.h>
#include <google/protobuf/text_format.h>
#include <algorithm>

namespace qrisp {

Structure::Structure(const Structure& rna) {
  size_ = rna.Size();
  sequence_.assign(rna.sequence_.cbegin(), rna.sequence_.cend());
  pairings_.assign(rna.pairings_.cbegin(), rna.pairings_.cend());
  quality_.assign(rna.quality_.cbegin(), rna.quality_.cend());
}

Structure::Structure(const string& brackets, const string& seq,
                     const vector<score_t>& qual) {
  size_ = brackets.size() + 1;
  LoadPairingsFromString(brackets);
  LoadSequenceFromString(seq);
  if (qual.size() > 0) {
    quality_.assign(qual.cbegin(), qual.cend());
    quality_.emplace(quality_.begin(), INVALID_QUALITY);
  }
}

bool Structure::InitializeFromAsciiPb(const string& input) {
  StructureMessage m;
  ::google::protobuf::TextFormat::ParseFromString(input, &m);
  return InitializeFromProto(m);
}

bool Structure::InitializeFromProto(const StructureMessage& rna) {
  size_ = rna.rows_size() + 1;
  if (size_ < MIN_RNA_SIZE || size_ > MAX_RNA_SIZE) {
    return false;
  }
  pairings_.assign(size_, IDX_NOT_SET);
  sequence_.assign(size_, INVALID_BASE);
  quality_.assign(size_, 0.0);
  for (const auto& row : rna.rows()) {
    pairings_[row.pos()] = row.pair();
    // if (row.pair() > != 0) {
    //  pairings_[row.pair()] = row.pos();
    //}
    sequence_[row.pos()] = row.base();
    quality_[row.pos()] = row.score();
  }
  return true;
}

bool Structure::LoadPairingsFromString(const string& input) {
  pairings_.assign(input.size() + 1, IDX_NOT_SET);
  BracketsToPairings(input, &pairings_);
  return true;
}

bool Structure::LoadPairingsFromVector(const vector<idx_t>& input) {
  CHECK(input.size() + 1 == size_);
  pairings_.assign(input.cbegin(), input.cend());
  pairings_.emplace(pairings_.begin(), IDX_NOT_SET);
  return true;
}

bool Structure::LoadSequenceFromString(const string& input) {
  if (input.size() + 1 == size_) {
    sequence_.assign(size_, IDX_NOT_SET);
    for (size_t i = 0; i < input.size(); i++) {
      idx_t base = input[i] - 65;
      CHECK(base <= 20);
      sequence_[i + 1] = char2int[base];
    }
    return true;
  }
  return false;
}

bool Structure::LoadSequenceFromVector(const vector<idx_t>& input) {
  CHECK(input.size() + 1 == size_);
  sequence_.assign(input.begin(), input.end());
  sequence_.emplace(sequence_.begin(), INVALID_BASE);
  return true;
}

void Structure::ConvertToProto(StructureMessage* rna) const {
  for (size_t i = 1; i < pairings_.size(); i++) {
    auto* row = rna->add_rows();
    row->set_pos(i);
    row->set_pair(pairings_[i]);
    row->set_base(sequence_[i]);
    row->set_score(quality_[i]);
  }
}

Tuple Structure::GetBasesAt(const Tuple& t) const {
  // const idx_t len= this->Size();
  idx_t i, j, k, l;
  std::tie(i, j, k, l) = t;
  // if (!(i > 0 && (i < len || i == IDX_NOT_SET) &&
  //      j > 0 && (j < len || j == IDX_NOT_SET) &&
  //      k > 0 && (k < len || k == IDX_NOT_SET) &&
  //      l > 0 && (l < len || l == IDX_NOT_SET))) {
  //  printf("Invalid indices: %d %d %d %d\n", i, j, k, l);
  //  exit(EXIT_FAILURE);
  //}
  Tuple ret(i == IDX_NOT_SET ? 100 : sequence_[i],
            j == IDX_NOT_SET ? 100 : sequence_[j],
            k == IDX_NOT_SET ? 100 : sequence_[k],
            l == IDX_NOT_SET ? 100 : sequence_[l]);
  return ret;
}

void Structure::GetPairings(vector<idx_t>* p) const {
  p->clear();
  p->insert(p->begin(), pairings_.cbegin(), pairings_.cend());
}

void Structure::SetPairings(const vector<idx_t>& p) {
  pairings_.clear();
  pairings_.insert(pairings_.begin(), p.cbegin(), p.cend());
}

void Structure::SetQuality(const vector<score_t>& q) {
  quality_.clear();
  quality_.insert(quality_.begin(), q.cbegin(), q.cend());
}

void Structure::GetQuality(vector<score_t>* q) const {
  q->clear();
  q->insert(q->begin(), quality_.cbegin(), quality_.cend());
}

Structure& Structure::operator=(const Structure& s) {
  size_ = s.size_;
  pairings_.clear();
  pairings_.insert(pairings_.begin(), s.pairings_.cbegin(), s.pairings_.cend());
  sequence_.clear();
  sequence_.insert(sequence_.begin(), s.sequence_.cbegin(), s.sequence_.cend());
  quality_.clear();
  quality_.insert(quality_.begin(), s.quality_.cbegin(), s.quality_.cend());
  return *this;
}

// This function checks whether we have any pseudo-knots.
bool Structure::ContainsPseudoKnot() const {
  vector<std::tuple<idx_t, idx_t>> all_pairs;
  for (size_t i = 1; i < this->Size(); i++) {
    if (pairings_[i] != 0 && i < pairings_[i]) {
      all_pairs.push_back(std::make_tuple(i, pairings_[i]));
    }
  }
  for (size_t k = 1; k < all_pairs.size(); k++) {
    idx_t i, j, ip, jp;
    std::tie(i, j) = all_pairs[k - 1];
    std::tie(ip, jp) = all_pairs[k];
    if ((i <= ip && ip <= j) || (ip <= j && j <= jp)) {
      return true;
    }
  }
  return false;
}

// This function performs some basic sanity checks.
bool Structure::IsValidStructure() const {
  const size_t len = this->Size();
  // Just following a convention.
  CHECK(this->pairings_[0] == IDX_NOT_SET);
  // At least one position has to be in the structure.
  if (len < MIN_RNA_SIZE || len > MAX_RNA_SIZE) {
    return false;
  }
  // Perform some sanity checks on the structure.
  for (size_t i = 1; i < len; i++) {
    const idx_t j = pairings_[i];
    const idx_t ip = pairings_[j];
    if (j >= len) {
      return false;
    }
    if (j != 0 && i != ip) {
      return false;
    }
  }
  return true;
}

// This function calulates the base pairs of all 0-,1- and 2-loops of a given
// structure.
void Structure::CalculateSubstructure(
    vector<Substructure>* result, vector<bool>* accessible_positions) const {
  result->clear();
  const size_t len = this->Size();
  accessible_positions->resize(len, false);
  if (!this->IsValidStructure()) {
    return;
  }
  for (size_t i = 1; i < len; i++) {
    const idx_t pos = i;
    const idx_t pair = pairings_[pos];
    if (pair != 0 && pos < pair) {
      Substructure sub;
      sub.SetOuterPair(pos, pair);
      (*accessible_positions)[pos] = true;
      (*accessible_positions)[pair] = true;
      VLOG(1) << "Outer pair: " << pos << ", " << pair << endl;
      for (size_t iter_pos = pos + 1; iter_pos < pair;) {
        idx_t iter_pos_pair = pairings_[iter_pos];
        (*accessible_positions)[iter_pos] = true;
        if (iter_pos >= iter_pos_pair) {
          iter_pos++;
          continue;
        }
        VLOG(1) << "\tinner pair: " << iter_pos << ", " << iter_pos_pair
                << endl;
        // We have no pairing at current position.
        if (iter_pos_pair == 0) {
          iter_pos++;
          continue;
        }
        // We have a pairing at the given position.
        if (iter_pos_pair != 0) {
          sub.AddInnerPair(iter_pos, iter_pos_pair);
          iter_pos = iter_pos_pair + 1;
          continue;
        }
      }
      result->push_back(sub);
    }
  }
}

double PairwiseLoss(const Structure& a, const Structure& b) {
  CHECK(a.Size() == b.Size());
  double cur_loss = 0.0;
  for (size_t i = 1; i < a.Size(); i++) {
    if (a(i) != b(i)) {
      cur_loss += 1.0;
    }
  }
  return cur_loss;
}

double HammingLoss(const Structure& a, const Structure& b) {
  CHECK(a.Size() == b.Size());
  double cur_loss = 0.0;
  for (size_t i = 1; i < a.Size(); i++) {
    if ((a(i) == 0 && b(i) != 0) || (a(i) != 0 && b(i) == 0)) {
      cur_loss += 1.0;
    }
  }
  return cur_loss;
}

void PairingsToBrackets(const vector<idx_t>& pairings, string* brackets) {
  brackets->clear();
  brackets->resize(pairings.size() - 1, '.');
  for (size_t i = 1; i < pairings.size(); i++) {
    if (i < pairings[i] && pairings[i] != 0) {
      (*brackets)[i - 1] = '(';
      (*brackets)[pairings[i] - 1] = ')';
    }
  }
}

void BracketsToPairings(const string& brackets, vector<idx_t>* pairings) {
  vector<idx_t> opening;
  for (size_t i = 0; i < brackets.size(); i++) {
    if (brackets[i] == '(') {
      opening.push_back(i + 1);
      continue;
    }
    if (brackets[i] == ')') {
      (*pairings)[opening.back()] = i + 1;
      (*pairings)[i + 1] = opening.back();
      opening.pop_back();
      continue;
    }
    (*pairings)[i + 1] = 0;
  }
}

bool operator==(const Structure& a, const Structure& b) {
  if (a.Size() != b.Size()) {
    return false;
  }
  for (size_t i = 1; i < a.Size(); i++) {
    if (a(i) != b(i) || a[i] != b[i] || (a.quality(i) != b.quality(i))) {
      return false;
    }
  }
  return true;
}

bool operator==(const Substructure& a, const Substructure& b) {
  if (a.GetOuterPair() == b.GetOuterPair() &&
      a.GetNumberOfInnerPairs() == b.GetNumberOfInnerPairs()) {
    for (size_t i = 0; i < a.GetNumberOfInnerPairs(); i++) {
      if (a.GetInnerPairAt(i) != b.GetInnerPairAt(i)) {
        return false;
      }
    }
    return true;
  }
  return false;
}

Substructure::Substructure(const IntPair& outer, const vector<IntPair>& inner) {
  outer_pair_ = outer;
  inner_pairs_.insert(inner_pairs_.begin(), inner.cbegin(), inner.cend());
}

}  // namespace qrisp
