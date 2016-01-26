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

#ifndef QRISP_UTILS_H__
#define QRISP_UTILS_H__

#include <glog/logging.h>
#include <cfloat>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace qrisp {

using namespace std;

typedef int32_t idx_t;
typedef double score_t;
typedef const idx_t& CIR;

class Structure;

typedef map<string, Structure> Dataset;

constexpr idx_t MIN_RNA_SIZE = 6;
constexpr idx_t MAX_RNA_SIZE = 10000;

// constexpr idx_t IDX_NOT_SET = 20000;
// constexpr idx_t INVALID_BASE = 30000;
// constexpr idx_t NDX = 40000;

constexpr score_t LOG_ZERO = -numeric_limits<score_t>::max();
constexpr idx_t IDX_NOT_SET = -1;
constexpr idx_t INVALID_BASE = -2;
constexpr idx_t NDX = -3;
constexpr score_t INVALID_QUALITY = -1.0;

typedef pair<idx_t, idx_t> IntPair;
typedef std::tuple<idx_t, idx_t, idx_t, idx_t> Tuple;
typedef std::unordered_map<idx_t, score_t> FeatureVec;
typedef vector<pair<idx_t, score_t> > IndexVec;
typedef vector<score_t> ScoreVec;

// const ScoreVec limits = { 0.0, 0.25, 0.5, 0.75, 1.1 };
const ScoreVec limits = {0.0, 0.25, 0.5, 1.5, 3.0};

// Define mapping from characters to indices.
static int char2int[21] = {0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0, 0, 3, 3};

constexpr int b0 = 1;
constexpr int b1 = 4;
constexpr int b2 = 16;
constexpr int b3 = 64;

// This defines a mapping from 4 tuples representing stacked base pairs to an
// integer representing the feature offset for that particular pair of base
// pairs.
static int mapping[] = {
    0,  1,  2,   3,   4,  5,  6,   7,   8,  9,  10,  11,  12, 13, 14,  15,
    4,  16, 17,  18,  19, 20, 21,  22,  23, 24, 25,  26,  27, 28, 29,  30,
    8,  31, 32,  33,  23, 34, 35,  36,  37, 38, 39,  40,  41, 42, 43,  44,
    12, 45, 46,  47,  27, 48, 49,  50,  41, 51, 52,  53,  54, 55, 56,  57,
    1,  58, 59,  60,  16, 61, 62,  63,  31, 64, 65,  66,  45, 67, 68,  69,
    5,  61, 70,  71,  20, 72, 73,  74,  34, 75, 76,  77,  48, 78, 79,  80,
    9,  64, 81,  82,  24, 75, 83,  84,  38, 85, 86,  87,  51, 88, 89,  90,
    13, 67, 91,  92,  28, 78, 93,  94,  42, 88, 95,  96,  55, 97, 98,  99,
    2,  59, 100, 101, 17, 70, 102, 103, 32, 81, 104, 105, 46, 91, 106, 107,
    6,  62, 102, 108, 21, 73, 109, 110, 35, 83, 111, 112, 49, 93, 113, 114,
    10, 65, 104, 115, 25, 76, 111, 116, 39, 86, 117, 118, 52, 95, 119, 120,
    14, 68, 106, 121, 29, 79, 113, 122, 43, 89, 119, 123, 56, 98, 124, 125,
    3,  60, 101, 126, 18, 71, 108, 127, 33, 82, 115, 128, 47, 92, 121, 129,
    7,  63, 103, 127, 22, 74, 110, 130, 36, 84, 116, 131, 50, 94, 122, 132,
    11, 66, 105, 128, 26, 77, 112, 131, 40, 87, 118, 133, 53, 96, 123, 134,
    15, 69, 107, 129, 30, 80, 114, 132, 44, 90, 120, 134, 57, 99, 125, 135};

inline idx_t CalculateOffset(const Tuple& t) {
  idx_t i, j, k, l;
  std::tie(i, j, k, l) = t;
  CHECK(i != 100);
  CHECK(j != 100);
  CHECK(k != 100);
  CHECK(l != 100);
  idx_t index1 = i * b3 + j * b2 + k * b1 + l * b0;
  idx_t index2 = l * b3 + k * b2 + j * b1 + i * b0;
  if (index1 <= index2) {
    return mapping[index1];
  }
  return mapping[index2];
}

inline idx_t CalculateFullOffset(const Tuple& t) {
  idx_t i, j, k, l;
  std::tie(i, j, k, l) = t;
  return i * b3 + j * b2 + k * b1 + l * b0;
}

inline idx_t CalculateLeftOffset(const Tuple& t) {
  int i, j, k;
  std::tie(i, j, k, std::ignore) = t;
  return i * b2 + j * b1 + k * b0;
}

inline idx_t CalculateRightOffset(const Tuple& t) {
  int i, j, l;
  std::tie(i, j, std::ignore, l) = t;
  return i * b2 + j * b1 + l * b0;
}

auto mt = [](idx_t i, idx_t j, idx_t k, idx_t l) {
  return std::make_tuple(i, j, k, l);
};

auto mt3 = [](idx_t i, idx_t j, idx_t k) {
  return std::make_tuple(i, j, k, IDX_NOT_SET);
};

auto mt2 = [](idx_t i, idx_t j) { return std::make_tuple(i, j, NDX, NDX); };

auto mt1 = [](idx_t i) {
  return std::make_tuple(i, IDX_NOT_SET, IDX_NOT_SET, IDX_NOT_SET);
};

auto mp = [](idx_t i, score_t s) -> pair<idx_t, score_t> {
  return std::make_pair(i, s);
};

// auto mpii = [] (idx_t i, idx_t j) -> pair<idx_t, idx_t> {
//  return make_pair(i, j);
//};

template <typename T>
class Array {
 public:
  Array(int s) : size_(s) { data_ = new T[size_]; }

  Array(int s, T init_elem) : size_(s) {
    data_ = new T[size_];
    Init(init_elem);
  }

  virtual ~Array() { delete[] data_; }

  void Init(T init_elem) {
    for (int i = 0; i < size_; i++) data_[i] = init_elem;
  }

  T& Get(const int& idx) const {
    CHECK(idx < size_);
    return data_[idx];
  }

  T& operator()(int idx) {
    CHECK(idx < size_);
    return data_[idx];
  }

  void SetElement(int idx, T elem) {
    CHECK(idx < size_);
    data_[idx] = elem;
  }

 private:
  int size_;
  T* data_;
};

template <typename T>
class Array3 : public Array<T> {
 public:
  Array3(int d1, int d2, int d3)
      : Array<T>(d1 * d2 * d3), d1_size(d1), d2_size(d2), d3_size(d3) {}

  Array3(int d1, int d2, int d3, T init_elem)
      : Array<T>(d1 * d2 * d3, init_elem),
        d1_size(d1),
        d2_size(d2),
        d3_size(d3) {}

  ~Array3() {}

  T& Get(const int& i, const int& j, const int& k) const {
    CHECK(i < d1_size && j < d2_size && k < d3_size);
    return Array<T>::Get(i + d1_size * (j + d2_size * k));
  }

  void Set(int i, int j, int k, T elem) {
    CHECK(i < d1_size && j < d2_size && k < d3_size);
    Array<T>::SetElement(i + d1_size * (j + d2_size * k), elem);
  }

  void add_to_element(int i, int j, int k, T elem) {
    CHECK(i < d1_size && j < d2_size && k < d3_size);
    T tmp = Array<T>::Get(i + d1_size * (j + d2_size * k));
    Array<T>::set_element(i + d1_size * (j + d2_size * k), tmp + elem);
  }

  void Init(T elem) { Array<T>::Init(elem); }

  T& operator()(const int& i, const int& j, const int& k) {
    if (i < d1_size && j < d2_size && k < d3_size) {
      return Array<T>::operator()(i + d1_size * (j + d2_size * k));
    }
    return Array<T>::operator()(0);
  }

  int GetSize1() const { return d1_size; }

  int GetSize2() const { return d2_size; }

  int GetSize3() const { return d3_size; }

 private:
  int d1_size;
  int d2_size;
  int d3_size;
};

typedef Array3<score_t> DPTable;
typedef std::tuple<double, double> Performance;

inline double F1Measure(const Performance& p) {
  double sens, spec;
  std::tie(sens, spec) = p;
  if (sens + spec > 0.0) {
    return 2 * (sens * spec) / (sens + spec);
  }
  return 0.0;
}

void SplitStringToInts(const string& input, vector<int>* output);

}  // namespace qrisp

#endif  // QRISP_UTILS_H__
