#include <limits.h>
#include <manyclaw/manyclaw.h>
#include "gtest/gtest.h"

// Test the indexer methods
TEST(DataStructuresTest, FieldIndexer) {
	FieldIndexer<10, 10, 2, 1> fi;

	EXPECT_EQ(196, fi.size());

	EXPECT_EQ(0, fi.idx(0,0));
	EXPECT_EQ(14, fi.idx(1,0));
	EXPECT_EQ(1, fi.idx(0,1));

	EXPECT_EQ(fi.idx(2,1), fi.up(1, 1));
	EXPECT_EQ(fi.idx(0,1), fi.down(1, 1));
	EXPECT_EQ(fi.idx(1,0), fi.left(1, 1));
	EXPECT_EQ(fi.idx(1,2), fi.right(1, 1));
}


