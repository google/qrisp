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

#include "rna-structure.h"

#include <glog/logging.h>
#include <algorithm>

namespace qrisp {

Structure::Structure(const Structure& rna) {
  size_ = rna.Size();
  has_quality_ = rna.has_quality_;
  CHECK(size_ >= MIN_RNA_SIZE);
  sequence_.assign(rna.sequence_.cbegin(), rna.sequence_.cend());
  pairings_.assign(rna.pairings_.cbegin(), rna.pairings_.cend());
  quality_.assign(rna.quality_.cbegin(), rna.quality_.cend());
}


bool Structure::InitializeFromProto(const StructureMessage& rna) {
  size_ = rna.length() + 1;
  has_quality_ = false;
  if (rna.base_size() > 0) {
    sequence_.assign(rna.base().begin(), rna.base().end());
  } else if (!rna.sequence().empty()) {
    const auto seq = rna.sequence();
    sequence_.assign(seq.size(), IDX_NOT_SET);
    for (int i = 0; i < seq.size(); i++) {
      idx_t base = seq[i] - 65;
      if (base > 20) {
        return false;
      }
      sequence_[i + 1] = char2int[base];
    }
  } else {
    return false;
  }
  // Initialize pairings with zero everywhere. Position zero is set to a dummy
  // value.
  pairings_.assign(rna.length(), 0);
  pairings_[0] = IDX_NOT_SET;
  if (rna.pair_size() > 0) {
    for (const auto& pair : rna.pair()) {
      pairings_.assign(rna.length(), 0);
      pairings_[pair.lo()] = pair.hi();
      pairings_[pair.hi()] = pair.lo();
    }
  } else if (!rna.brackets().empty()) {
    BracketsToPairings(rna.brackets(), &pairings_);
  } else {
    return false;
  }
  if (rna.confidence_size() + 1 == size_) {

  }
  pairings_.assign(rna.length(), 0);
  quality_.assign(rna.confidence().begin(), rna.confidence().end());
  if (rna.confidence_size() == pairings_.size()) {
    has_quality_ = true;
  }
  return true;
}

bool Structure::LoadPairingsFromString(const string& input) {
  pairings_.assign(input.size() + 1, IDX_NOT_SET);
  BracketsToPairings(input, &pairings_);
  return true;
}

bool Structure::LoadPairingsFromVector(const vector<idx_t>& input) {
  pairings_.assign(input.cbegin(), input.cend());
  return true;
}

bool Structure::LoadSequenceFromString(const string& input) {
  size_ = input.size() + 1;
  CHECK(size_ >= MIN_RNA_SIZE);
  CHECK(size_ < MAX_RNA_SIZE);
  sequence_.assign(size_, IDX_NOT_SET);
  for (int i = 0; i < input.size(); i++) {
    idx_t base = input[i] - 65;
    CHECK(base <= 20);
    sequence_[i + 1] = char2int[base];
  }
  return true;
}

bool Structure::LoadSequenceFromVector(const vector<idx_t>& input) {
  sequence_.assign(input.begin(), input.end());
  return true;
}

bool Structure::Initialize(const string& brackets, const string& seq,
                           const vector<score_t>& qual) {
  size_ = brackets.size() + 1;
  LoadPairingsFromString(brackets);
  LoadSequenceFromString(seq);
  if (qual.size() > 0) {
    quality_.assign(qual.cbegin(), qual.cend());
    has_quality_ = true;
  }
  if (pairings_.size() == size_ && sequence_.size() == size_ &&
      (quality_.size() == size_ || quality_.size() == 0)) {
    return true;
  }
  return false;
}

//Structure::Structure(const vector<idx_t>& pairings, const string& seq,
//                     const vector<score_t>& qual)
//  pairings_.assign(pairings.cbegin(), pairings.cend());
//  if (!qual.empty()) {
//    if (size_ != qual.size() + 1) {
//      LOG(ERROR) << "Quality vector has wrong size!";
//      exit(EXIT_FAILURE);
//    }
//    quality_.assign(qual.cbegin(), qual.cend());
//    quality_.emplace(quality_.begin(), INVALID_QUALITY);
//    has_quality_ = true;
//  }
//}


bool Structure::ConvertToProto(StructureMessage* rna) const {
  rna->set_length(pairings_.size() - 1);
  for (int i=1; i < pairings_.size(); i++) {
    if (pairings_[i] != 0 && i < pairings_[i]) {
      auto p = rna->add_pair();
      p->set_lo(i);
      p->set_lo(pairings_[i]);
    }
  }
  rna->mutable_base()->Resize(size_, -1);
  std::copy(sequence_.begin(), sequence_.end(), rna->mutable_base()->begin());
  rna->mutable_confidence()->Resize(size_, -1);
  std::copy(quality_.begin(), quality_.end(),
            rna->mutable_confidence()->begin());
  return true;
}

Tuple Structure::GetBasesAt(const Tuple& t) const {
  // const idx_t len= this->Size();
  int i, j, k, l;
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
  has_quality_ = true;
  quality_.clear();
  quality_.insert(quality_.begin(), q.cbegin(), q.cend());
}

void Structure::GetQuality(vector<score_t>* q) const {
  q->clear();
  q->insert(q->begin(), quality_.cbegin(), quality_.cend());
}

Structure& Structure::operator=(const Structure& s) {
  size_ = s.size_;
  has_quality_ = s.has_quality_;
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
  const int len = this->Size();
  // Just following a convention.
  CHECK(this->pairings_[0] == IDX_NOT_SET);
  // At least one position has to be in the structure.
  CHECK(len >= MIN_RNA_SIZE);
  CHECK(len < MAX_RNA_SIZE);
  // Perform some sanity checks on the structure.
  for (int i = 1; i < len; i++) {
    if (pairings_[i] >= len) {
      return false;
    }
  }
  // return !this->ContainsPseudoKnot();
  return true;
}


// This function calulates the base pairs of all 0-,1- and 2-loops of a given
// structure.
void Structure::CalculateSubstructure(
    vector<Substructure>* result, vector<bool>* accessible_positions) const {
  result->clear();
  const int len = this->Size();
  accessible_positions->resize(len, false);
  if (!this->IsValidStructure()) {
    return;
  }
  for (int i = 1; i < len; i++) {
    const idx_t pos = i;
    const idx_t pair = pairings_[pos];
    if (pair != 0 && pos < pair) {
      Substructure sub;
      sub.SetOuterPair(pos, pair);
      (*accessible_positions)[pos] = true;
      (*accessible_positions)[pair] = true;
      VLOG(1) << "Outer pair: " << pos << ", " << pair << endl;
      for (int iter_pos = pos + 1; iter_pos < pair;) {
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
  for (int i = 1; i < a.Size(); i++) {
    if (a(i) != b(i)) {
      cur_loss += 1.0;
    }
  }
  return cur_loss;
}

double HammingLoss(const Structure& a, const Structure& b) {
  CHECK(a.Size() == b.Size());
  double cur_loss = 0.0;
  for (int i = 1; i < a.Size(); i++) {
    if ((a(i) == 0 && b(i) != 0) || (a(i) != 0 && b(i) == 0)) {
      cur_loss += 1.0;
    }
  }
  return cur_loss;
}

void PairingsToBrackets(const vector<idx_t>& pairings, string* brackets) {
  brackets->clear();
  brackets->resize(pairings.size() - 1, '.');
  for (int i = 1; i < pairings.size(); i++) {
    if (i < pairings[i] && pairings[i] != 0) {
      (*brackets)[i - 1] = '(';
      (*brackets)[pairings[i] - 1] = ')';
    }
  }
}

void BracketsToPairings(const string& brackets, vector<idx_t>* pairings) {
  vector<idx_t> opening;
  for (idx_t i = 0; i < brackets.size(); i++) {
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
  for (int i = 1; i < a.Size(); i++) {
    if (a(i) != b(i) || a[i] != b[i] ||
        (a.HasQuality() && b.HasQuality() && a.quality(i) != b.quality(i))) {
      return false;
    }
  }
  return true;
}

bool operator==(const Substructure& a, const Substructure& b) {
  if (a.GetOuterPair() == b.GetOuterPair() &&
      a.GetNumberOfInnerPairs() == b.GetNumberOfInnerPairs()) {
    for (int i = 0; i < a.GetNumberOfInnerPairs(); i++) {
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
