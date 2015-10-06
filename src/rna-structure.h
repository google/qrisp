#ifndef QRISP_RNA_STRUCTURE_H__
#define QRISP_RNA_STRUCTURE_H__

#include <glog/logging.h>
#include "proto/structure.pb.h"
#include "utils.h"

#include <string>
#include <vector>

namespace qrisp {

using namespace std;

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

// This class implements RNA secondary structures.
//
// AAGCTATGCCA sequence
// ..((...)).. structure
// q1q2q3q4... quality
//
// The size of the structure is defined as the size of the sequence + 1.
class Structure {
 public:
  // This constructor expects a RNA structure in bracket notation form.
  Structure(const string& brackets, const string& seq,
            const vector<score_t>& qual);

  // The following constructors expects an RNA structure given as a vector of
  // pairs.
  Structure(const vector<idx_t>& pairings, const string& seq,
            const vector<score_t>& qual);

  // Constructor from RNASStructure protocol buffer.
  Structure(const RNASStructure& rna);

  // Copy constructor for Structure objects themselves.
  Structure(const Structure& rna);

  // Convert the structure back to the proto representation.
  void ConvertToProto(RNASStructure* rna) const;

  // Field / element accessing methods.
  inline idx_t GetSize() const { return size_; }

  // Returns pair information for position i.
  inline idx_t operator()(idx_t i) const {
    CHECK(0 < i && i < pairings_.size());
    return pairings_[i];
  }

  // Returns sequence information for position i.
  idx_t operator[](idx_t i) const;

  // Returns quality information for position i.
  score_t quality(idx_t i) const;

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

  inline bool HasQuality() const { return has_quality_; }

  // Performs some sanity checks.
  bool ContainsPseudoKnot() const;

  bool IsValidStructure() const;

  // Compare two structures.
  friend bool operator==(const Structure& a, const Structure& b);

 private:
  void AuxConstructor(const string& seq);

  vector<idx_t> pairings_;
  vector<idx_t> sequence_;
  vector<score_t> quality_;

  idx_t size_;

  bool has_quality_;

  // DISALLOW_IMPLICIT_CONSTRUCTORS(Structure);
};

inline idx_t Structure::operator[](idx_t i) const {
  CHECK(0 < i && i < sequence_.size());
  return sequence_[i];
}

inline double Structure::quality(idx_t i) const {
  CHECK(0 < i && i < quality_.size());
  return quality_[i];
}

void PairingsToBrackets(const vector<idx_t>& pairings, string* brackets);

void BracketsToPairings(const string& brackets, vector<idx_t>* pairings);

double PairwiseLoss(const Structure& a, const Structure& b);

double HammingLoss(const Structure& a, const Structure& b);

}  // namespace qrisp
#endif  // QRISP_RNA_STRUCTURE_H__
