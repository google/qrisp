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
//
// This library contains all classes and algorithms needed for the
// representation of one or more RNA secondary structures.

#ifndef QRISP_STRUCTURE_H__
#define QRISP_STRUCTURE_H__

#include <glog/logging.h>
#include <string>
#include <vector>

#include "qrisp/proto/structure.pb.h"
#include "utils.h"

namespace qrisp {

// Forward declaration. Definition below.
class Substructure;

// This class implements one RNA secondary structure.
//
// AAGCTATGCCA sequence
// ..((...)).. structure
// q1q2q3q4... quality
//
//
// The size of the structure is defined as the size of the sequence + 1.
class Structure {
 public:
  // Initialize an emtpy structure.
  Structure() : size_(0) {}

  // Copy constructor for Structure objects themselves.
  Structure(const Structure& rna);

  // Initialize the structure directly from string representations of the data.
  Structure(const string& brackets, const string& seq,
            const vector<score_t>& qual);

  // The size of the structure is defined as the size of the sequence + 1.
  inline idx_t Size() const { return size_; }

  // Operator () returns pair information for position i.
  inline idx_t operator()(idx_t i) const {
    CHECK(0 < i && i < pairings_.size());
    return pairings_[i];
  }

  // Operator [] returns nucleotide information for position i.
  inline idx_t operator[](idx_t i) const {
    CHECK(0 < i && i < sequence_.size());
    return sequence_[i];
  }

  inline double quality(idx_t i) const {
    LOG(INFO) << "quality idx: " << i;
    CHECK(0 < i && i < quality_.size());
    return quality_[i];
  }

  //
  Tuple GetBasesAt(const Tuple& t) const;

  //
  void GetPairings(vector<idx_t>* p) const;

  void SetPairings(const vector<idx_t>& p);
  //
  void GetQuality(vector<score_t>* q) const;

  void SetQuality(const vector<score_t>& q);

  //
  Structure& operator=(const Structure&);

  // Calculate some features of the current structure.
  void CalculateSubstructure(vector<Substructure>* result,
                             vector<bool>* accessible_positions) const;

  // Performs some sanity checks.
  bool ContainsPseudoKnot() const;

  bool IsValidStructure() const;

  // Compare two structures.
  friend bool operator==(const Structure& a, const Structure& b);

  // Initialize the Structure object with information coming from the proto.
  bool InitializeFromProto(const StructureMessage& rna);

  bool InitializeFromAsciiPb(const string& input);

  // Convert the structure back to the proto representation.
  void ConvertToProto(StructureMessage* rna) const;

 private:
  bool LoadPairingsFromString(const string& input);
  bool LoadPairingsFromVector(const vector<idx_t>& input);
  bool LoadSequenceFromString(const string& input);
  bool LoadSequenceFromVector(const vector<idx_t>& input);

  vector<idx_t> pairings_;
  vector<idx_t> sequence_;
  vector<score_t> quality_;

  idx_t size_;
};

// This constructor expects a RNA structure in bracket notation form.
// Structure(const string& brackets, const string& seq, const vector<score_t>&
// qual);

// The following constructors expects an RNA structure given as a vector of
// pairs.
//  Structure(const vector<idx_t>& pairings, const string& seq,
//           const vector<score_t>& qual);

void PairingsToBrackets(const vector<idx_t>& pairings, string* brackets);

void BracketsToPairings(const string& brackets, vector<idx_t>* pairings);

// In the following we define two loss functions:
// * Pairwise loss, and
// * Hamming loss.
//
// The distinction between pairwise and Hamming loss can be explained by the
// following example. Let a and b be defined as follows:
// a := .((....))
// b := (.(....))
// Then the pairwise loss would be one while the Hamming loss would be two.

// Calculate the pairwise loss between structures a and b.
double PairwiseLoss(const Structure& a, const Structure& b);

// Calculate the Hamming loss between structures a and b.
double HammingLoss(const Structure& a, const Structure& b);

// This class is used to store information about substructures occurring in an
// RNA secondary structure.
class Substructure {
 public:
  Substructure() : outer_pair_(make_pair(IDX_NOT_SET, IDX_NOT_SET)) {}

  Substructure(const IntPair& outer, const vector<IntPair>& inner);

  ~Substructure() {}

  IntPair GetOuterPair() const { return outer_pair_; }

  void SetOuterPair(const idx_t& i, const idx_t& j) {
    outer_pair_.first = i;
    outer_pair_.second = j;
  }

  void AddInnerPair(const idx_t& i, const idx_t& j) {
    inner_pairs_.push_back(make_pair(i, j));
  }

  idx_t GetNumberOfInnerPairs() const { return inner_pairs_.size(); }

  IntPair GetInnerPairAt(const idx_t& pos) const {
    assert(pos < inner_pairs_.size());
    return inner_pairs_[pos];
  }

 private:
  IntPair outer_pair_;

  vector<IntPair> inner_pairs_;

  friend bool operator==(const Substructure& a, const Substructure& b);
};

}  // namespace qrisp
#endif  // QRISP_STRUCTURE_H__
