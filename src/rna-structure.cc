#include "rna-structure.h"

#include <glog/logging.h>
#include <algorithm>

namespace qrisp {

void Structure::AuxConstructor(const string& seq) {
  size_ = seq.size() + 1;
  CHECK(size_ >= MIN_RNA_SIZE);
  CHECK(size_ < MAX_RNA_SIZE);
  sequence_.assign(size_, IDX_NOT_SET);
  for (int i = 0; i < seq.size(); i++) {
    idx_t base = seq[i] - 65;
    CHECK(base <= 20);
    sequence_[i + 1] = char2int[base];
  }
  pairings_.assign(size_, IDX_NOT_SET);
}

Structure::Structure(const string& brackets, const string& seq,
                     const vector<score_t>& qual)
    : has_quality_(false) {
  AuxConstructor(seq);
  CHECK(size_ == brackets.size() + 1);
  BracketsToPairings(brackets, &pairings_);
  if (!qual.empty()) {
    if (size_ != qual.size() + 1) {
      LOG(ERROR) << "Quality vector has wrong size!";
      exit(EXIT_FAILURE);
    }
    quality_.assign(qual.cbegin(), qual.cend());
    quality_.emplace(quality_.begin(), INVALID_QUALITY);
    has_quality_ = true;
  }
}

Structure::Structure(const vector<idx_t>& pairings, const string& seq,
                     const vector<score_t>& qual)
    : has_quality_(false) {
  AuxConstructor(seq);
  pairings_.assign(pairings.cbegin(), pairings.cend());
  if (!qual.empty()) {
    if (size_ != qual.size() + 1) {
      LOG(ERROR) << "Quality vector has wrong size!";
      exit(EXIT_FAILURE);
    }
    quality_.assign(qual.cbegin(), qual.cend());
    quality_.emplace(quality_.begin(), INVALID_QUALITY);
    has_quality_ = true;
  }
}

Structure::Structure(const RNASStructure& rna) {
  size_ = rna.base_size();
  has_quality_ = false;
  sequence_.assign(rna.base().begin(), rna.base().end());
  pairings_.assign(rna.pair().begin(), rna.pair().end());
  quality_.assign(rna.confidence().begin(), rna.confidence().end());
  if (rna.confidence_size() == pairings_.size()) {
    has_quality_ = true;
  }
}

Structure::Structure(const Structure& rna) {
  size_ = rna.GetSize();
  has_quality_ = rna.has_quality_;
  CHECK(size_ >= MIN_RNA_SIZE);
  sequence_.assign(rna.sequence_.cbegin(), rna.sequence_.cend());
  pairings_.assign(rna.pairings_.cbegin(), rna.pairings_.cend());
  quality_.assign(rna.quality_.cbegin(), rna.quality_.cend());
}

void Structure::ConvertToProto(RNASStructure* rna) const {
  rna->mutable_base()->Resize(size_, -1);
  std::copy(sequence_.begin(), sequence_.end(), rna->mutable_base()->begin());
  rna->mutable_pair()->Resize(size_, -1);
  std::copy(pairings_.begin(), pairings_.end(), rna->mutable_pair()->begin());
  rna->mutable_confidence()->Resize(size_, -1);
  std::copy(quality_.begin(), quality_.end(),
            rna->mutable_confidence()->begin());
}

Tuple Structure::GetBasesAt(const Tuple& t) const {
  // const idx_t len= this->GetSize();
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
  for (int i = 1; i < this->GetSize(); i++) {
    if (pairings_[i] != 0 && i < pairings_[i]) {
      all_pairs.push_back(std::make_tuple(i, pairings_[i]));
    }
  }
  for (int k = 1; k < all_pairs.size(); k++) {
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
  const int len = this->GetSize();
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
  const int len = this->GetSize();
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
  CHECK(a.GetSize() == b.GetSize());
  double cur_loss = 0.0;
  for (int i = 1; i < a.GetSize(); i++) {
    if (a(i) != b(i)) {
      cur_loss += 1.0;
    }
  }
  return cur_loss;
}

double HammingLoss(const Structure& a, const Structure& b) {
  CHECK(a.GetSize() == b.GetSize());
  double cur_loss = 0.0;
  for (int i = 1; i < a.GetSize(); i++) {
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
  if (a.GetSize() != b.GetSize()) {
    return false;
  }
  for (int i = 1; i < a.GetSize(); i++) {
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
