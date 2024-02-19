//
// Created by Oscar Tse on 31/12/2023.
//

#include <gtest/gtest.h>
#include "library.h"
#include <vector>
#include <string>
#include <cmath>

#include <fstream>

std::vector<Segment<uint64_t, double>> getFromRawString(const std::string &rawString) {
    std::vector<Segment<uint64_t, double>> segments;

    // Create a string stream from the raw string
    std::istringstream ss(rawString);

    // Read each line of the raw string
    std::string line;
    while (std::getline(ss, line)) {
        // Parse the line and create a Segment object
        std::istringstream line_ss(line);
        std::string x_start_str, slope_str, y_str;
        std::getline(line_ss, x_start_str, ',');
        std::getline(line_ss, slope_str, ',');
        std::getline(line_ss, y_str, ',');
        uint64_t x_start = std::stoull(x_start_str);
        double slope = std::stod(slope_str);
        double y = std::stod(y_str);
        Segment<uint64_t, double> segment(x_start, slope, y);

        // Add the Segment object to the vector
        segments.push_back(segment);
    }

    return segments;
}

TEST(PointTest, PointRetrival) {
    auto s = Point<double>(0.5, 0.5);
    EXPECT_DOUBLE_EQ(0.5, s.x);
    EXPECT_DOUBLE_EQ(0.5, s.y);
}

TEST(PointTest, DefaultPoint) {
    auto s = Point<double>();
    EXPECT_DOUBLE_EQ(0, s.x);
    EXPECT_DOUBLE_EQ(0, s.y);
}

TEST(PointTest, GetUpperLowerBound) {
    auto s = Point<double>(1, 1);
    EXPECT_DOUBLE_EQ(2, s.getUpperBound(1).y);
    EXPECT_DOUBLE_EQ(0, s.getLowerBound(1).y);
}


TEST(LineTest, DefaultLineTest) {
    auto s = Line<double>();
    EXPECT_DOUBLE_EQ(s.a1, 0);
    EXPECT_DOUBLE_EQ(s.a2, 0);
}

TEST(LineTest, LineConstructor) {
    auto s = Line<double>(1, 0);
    EXPECT_DOUBLE_EQ(1, s.a1);
    EXPECT_DOUBLE_EQ(0, s.a2);
}

TEST(LineTest, LineAbove) {
    auto s = Line<double>(1, 0);
    auto p1 = Point<double>(2, 3);
    EXPECT_EQ(s.above(p1), true);

    // same coord shd be false
    auto p3 = Point<double>(0, 0);
    EXPECT_EQ(s.above(p3), false);
}

TEST(LineTest, LineBelow) {
    auto s = Line<double>(1, 2);
    auto p1 = Point<double>(1, 1);
    EXPECT_EQ(s.below(p1), true);

    // same coord shd return false
    auto p2 = Point<double>(1, 3);
    EXPECT_EQ(s.below(p2), false);
}


TEST(LineTest, LineIntersection) {
    auto l1 = Line<double>(1, 0);
    auto l2 = Line<double>(-1, 0);
    auto i_p = l1.getIntersection(l2);
    EXPECT_DOUBLE_EQ(i_p.x, 0);
    EXPECT_DOUBLE_EQ(i_p.y, 0);

    auto l3 = Line<double>(1, 3);
    auto l4 = Line<double>(-1, -3);
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

TEST(PLRDataRepTest, TestEncodeDecode) {
    auto a1 = Segment<uint64_t, double>(1, 3.14, 5.67);
    auto a2 = Segment<uint64_t, double>(3, 6.78, 9);
    auto a3 = Segment<uint64_t, double>(8, 4.3443, 9.314);
    auto a4 = Segment<uint64_t, double>(18374830, 8.431413, 4332.124344);
    PLRDataRep<uint64_t, double> plrDataRep = PLRDataRep<uint64_t, double>(0.005);
    plrDataRep.Add(a1);
    plrDataRep.Add(a2);
    plrDataRep.Add(a3);
    plrDataRep.Add(a4);
    auto encodedStr = plrDataRep.Encode();
//    std::cout << encodedStr << std::endl;
    auto decodedObj = PLRDataRep<uint64_t, double>(encodedStr);
    EXPECT_DOUBLE_EQ(decodedObj.GetGamma(), 0.005);

    EXPECT_TRUE(decodedObj.GetSegs()[0] == a1);
    EXPECT_TRUE(decodedObj.GetSegs()[1] == a2);
    EXPECT_TRUE(decodedObj.GetSegs()[2] == a3);
    EXPECT_TRUE(decodedObj.GetSegs()[3] == a4);
}

TEST(PLRDataRepTest, TestBinarySearch) {
    auto str = "1,0.00205553,-0.00205553\n884,0.00217752,0.0750024\n1840,0.00155507,1.13784\n2957,0.00133835,2.04209\n4156,0.00272838,-3.34017\n5153,0.00315824,-6.27621\n6151,0.00184995,0.620625\n7103,0.00190864,0.441149\n8001,0.00249316,-3.9502\n8853,0.00232765,-2.60776\n9567,0.00175636,3.19647\n10741,0.00256332,-5.53276\n11542,0,24";
    auto res = getFromRawString(str);
    auto plrDataRep = PLRDataRep<uint64_t, double>(0.0005, res);
    auto test1 = plrDataRep.GetValue(6152);
    auto test2 = plrDataRep.GetValue(9661);
    auto test3 = plrDataRep.GetValue(1990);
    EXPECT_DOUBLE_EQ(test1.first, floor(0.00184995 * 6152 + 0.620625-plrDataRep.GetGamma()));
    EXPECT_DOUBLE_EQ(test2.first, floor(9661 * 0.00175636 + 3.19647-plrDataRep.GetGamma()));
    EXPECT_DOUBLE_EQ(test3.first, floor(0.00155507* 1990 + 1.13784-plrDataRep.GetGamma()));
}

TEST(PLRDataRepTest, testNull) {
    std::string str = "a";
    auto res = to_type<uint64_t>(str);
    EXPECT_EQ(res, 97);
}

TEST(StrToUint, testNull) {
    auto str = "";
    auto i = stringToNumber<uint64_t>(str);
    EXPECT_EQ(i, 0);
}

TEST(StrToUint, testAcharacter) {
    std::string str = "a";
//    std::cout <<str.length() << std::endl;
    auto i = stringToNumber<uint64_t>(str);
    EXPECT_EQ(i, 6989586621679009792);
}

TEST(StrToUint, testAllCharacter) {
    std::string str = "!!!!!!!!";
//    std::cout <<str.length() << std::endl;
    auto i = stringToNumber<uint64_t>(str);
    EXPECT_EQ(i, 2387225703656530209);
}

int main() {
    ::testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}