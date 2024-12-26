#include <gtest/gtest.h>
#include <BasisSet/slater_primitive.hpp>

TEST(SlaterPrimitive, Constructor) {
    SlaterPrimitive sp(1, 0, 0, 1.0);
    EXPECT_EQ(sp.n(), 1);
    EXPECT_EQ(sp.l(), 0);
    EXPECT_EQ(sp.m(), 0);
    EXPECT_EQ(sp.alpha(), 1.0);
    EXPECT_NEAR(sp.normalization(), 2.0, 1e-6);

    SlaterPrimitive sp1(2, 1, 1, 3.0);
    EXPECT_NEAR(sp1.normalization(), 18.0, 1e-6);

    SlaterPrimitive sp2(3, 2, 1, 3.0);
    EXPECT_NEAR(sp2.normalization(), 19.718012070, 1e-6);
}