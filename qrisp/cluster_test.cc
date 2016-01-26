#include "experimental/users/fdb/qrsp/cluster.h"
#include "experimental/users/fdb/qrsp/structure.h"

#include "testing/base/public/gunit.h"

namespace fdb {
namespace qrsp {
namespace {

class ClusterTest : public ::testing::Test {
 protected:
  ClusterTest() {}

  ~ClusterTest() override {}
};

TEST_F(ClusterTest, InputMultisetTooSmall) {
  vector<vector<int>> structures;
  vector<int> centroid;
  EXPECT_FALSE(CalculateCentroid(structures, &centroid));
}

TEST_F(ClusterTest, InputMultisetOnlyOneElement) {
  vector<vector<int>> structures(1, vector<int>(10, 0));
  const string brac = "(((...)))";
  BracketsToPairings(brac, &structures[0]);
  vector<int> centroid(structures[0].size(), 0);
  EXPECT_TRUE(CalculateCentroid(structures, &centroid));
  string centroid_brac;
  PairingsToBrackets(centroid, &centroid_brac);
  EXPECT_EQ(centroid_brac, brac);
}

TEST_F(ClusterTest, InputMultisetTwoElements) {
  vector<vector<int>> structures(2, vector<int>(10, 0));
  const string brac = "(((...)))";
  const string brac2 = "(.(...).)";
  BracketsToPairings(brac, &structures[0]);
  BracketsToPairings(brac2, &structures[1]);
  vector<int> centroid(structures[0].size(), 0);
  EXPECT_TRUE(CalculateCentroid(structures, &centroid));
  string centroid_brac;
  PairingsToBrackets(centroid, &centroid_brac);
  EXPECT_EQ(centroid_brac, brac2);
}

TEST_F(ClusterTest, SameInputSameOutput) {
  vector<vector<int>> structures(3, vector<int>(10, 0));
  const string brac = "(((...)))";
  BracketsToPairings(brac, &structures[0]);
  BracketsToPairings(brac, &structures[1]);
  BracketsToPairings(brac, &structures[2]);
  vector<int> centroid(structures[0].size(), 0);
  EXPECT_TRUE(CalculateCentroid(structures, &centroid));
  string centroid_brac;
  PairingsToBrackets(centroid, &centroid_brac);
  EXPECT_EQ(centroid_brac, brac);
}

TEST_F(ClusterTest, TieInOnePairing) {
  vector<vector<int>> structures(4, vector<int>(10, 0));
  BracketsToPairings("(((...)))", &structures[0]);
  BracketsToPairings("(.(...).)", &structures[1]);
  BracketsToPairings("(((...)))", &structures[2]);
  BracketsToPairings("(.(...).)", &structures[3]);
  vector<int> centroid(structures[0].size(), 0);
  EXPECT_TRUE(CalculateCentroid(structures, &centroid));
  string centroid_brac;
  PairingsToBrackets(centroid, &centroid_brac);
  EXPECT_EQ("(.(...).)", centroid_brac);
}

TEST_F(ClusterTest, TieInAllPairingsEmptyStructure) {
  vector<vector<int>> structures(4, vector<int>(10, 0));
  BracketsToPairings("((.....))", &structures[0]);
  BracketsToPairings("(.......)", &structures[1]);
  BracketsToPairings(".((...)).", &structures[2]);
  BracketsToPairings("..(...)..", &structures[3]);
  vector<int> centroid(structures[0].size(), 0);
  EXPECT_TRUE(CalculateCentroid(structures, &centroid));
  string centroid_brac;
  PairingsToBrackets(centroid, &centroid_brac);
  EXPECT_EQ(".........", centroid_brac);
}

TEST_F(ClusterTest, BestOfThree) {
  vector<vector<int>> structures(3, vector<int>(10, 0));
  const string brac = "(((...)))";
  BracketsToPairings(brac, &structures[0]);
  BracketsToPairings("(.(...).)", &structures[1]);
  BracketsToPairings(brac, &structures[2]);
  vector<int> centroid(structures[0].size(), 0);
  EXPECT_TRUE(CalculateCentroid(structures, &centroid));
  string centroid_brac;
  PairingsToBrackets(centroid, &centroid_brac);
  EXPECT_EQ(centroid_brac, brac);
}

TEST_F(ClusterTest, BestOfThreeAnotherOne) {
  vector<vector<int>> structures(3, vector<int>(10, 0));
  BracketsToPairings("(((...)))", &structures[0]);
  BracketsToPairings("(.(...).)", &structures[1]);
  BracketsToPairings("(.(...).)", &structures[2]);
  vector<int> centroid(structures[0].size(), 0);
  EXPECT_TRUE(CalculateCentroid(structures, &centroid));
  string centroid_brac;
  PairingsToBrackets(centroid, &centroid_brac);
  EXPECT_EQ(centroid_brac, "(.(...).)");
}

TEST_F(ClusterTest, BestOfSeven) {
  vector<vector<int>> structures(3, vector<int>(10, 0));
  const string brac = "(((...)))";
  BracketsToPairings(brac, &structures[0]);
  BracketsToPairings("(.(...).)", &structures[1]);
  BracketsToPairings(brac, &structures[2]);
  vector<int> centroid(structures[0].size(), 0);
  EXPECT_TRUE(CalculateCentroid(structures, &centroid));
  string centroid_brac;
  PairingsToBrackets(centroid, &centroid_brac);
  EXPECT_EQ(centroid_brac, brac);
}

}  // namespace
}  // namespace qrsp
}  // namespace fdb
