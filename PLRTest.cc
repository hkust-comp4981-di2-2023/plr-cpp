//
// Created by Oscar Tse on 31/12/2023.
//

#include <gtest/gtest.h>
#include "library.h"

TEST(PointTest, PointRetrival) {
    auto s = Point<double>(0.5,0.5);
    EXPECT_DOUBLE_EQ(0.5,s.x);
    EXPECT_DOUBLE_EQ(0.5,s.y);
}

TEST(PointTest, DefaultPoint) {
    auto s = Point<double>();
    EXPECT_DOUBLE_EQ(0,s.x);
    EXPECT_DOUBLE_EQ(0,s.y);
}

TEST(PointTest, GetUpperLowerBound) {
    auto s = Point<double>(1,1);
    EXPECT_DOUBLE_EQ(2,s.getUpperBound(1).y);
    EXPECT_DOUBLE_EQ(0, s.getLowerBound(1).y);
}


TEST(LineTest, DefaultLineTest) {
    auto s = Line<double>();
    EXPECT_DOUBLE_EQ(s.a1,0);
    EXPECT_DOUBLE_EQ(s.a2,0);
}

TEST(LineTest, LineConstructor) {
    auto s = Line<double>(1,0);
    EXPECT_DOUBLE_EQ(1,s.a1);
    EXPECT_DOUBLE_EQ(0,s.a2);
}

TEST(LineTest, LineAbove) {
    auto s = Line<double>(1,0);
    auto p1 = Point<double>(2,3);
    EXPECT_EQ(s.above(p1), true);

    // same coord shd be false
    auto p3 = Point<double>(0,0);
    EXPECT_EQ(s.above(p3), false);
}

TEST(LineTest, LineBelow) {
    auto s = Line<double>(1,2);
    auto p1 = Point<double>(1,1);
    EXPECT_EQ(s.below(p1), true);

    // same coord shd return false
    auto p2 = Point<double>(1,3);
    EXPECT_EQ(s.below(p2), false);
}


TEST(LineTest, LineIntersection) {
    auto l1 = Line<double>(1,0);
    auto l2 = Line<double>(-1,0);
    auto i_p = l1.getIntersection(l2);
    EXPECT_DOUBLE_EQ(i_p.x,0);
    EXPECT_DOUBLE_EQ(i_p.y,0);

    auto l3 = Line<double>(1,3);
    auto l4 = Line<double>(-1,-3);
    auto i_p2 = l3.getIntersection(l4);
    EXPECT_DOUBLE_EQ(i_p2.y, 0);
    EXPECT_DOUBLE_EQ(i_p2.x, -3);
}

TEST(LineTest, ThrowAssertWhenSlopeSame) {
    auto l1 = Line<double>(1, 0);
    auto l2 = Line<double>{l1};
    // Will die because of the same slope error
    ASSERT_DEATH(l1.getIntersection(l2), "");
}

int main() {
    ::testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}