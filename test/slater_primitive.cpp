#include <gtest/gtest.h>
#include <BasisSet/slater_primitive.hpp>
#include <stdexcept>

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

    EXPECT_THROW(
        {
            SlaterPrimitive sp3(0, 0, 0, 1.0);
        },
        std::invalid_argument
    );

    EXPECT_THROW(
        {
            SlaterPrimitive sp3(1, 2, 0, 1.0);
        },
        std::invalid_argument
    );

    EXPECT_THROW(
        {
            SlaterPrimitive sp3(1, 0, 3, 1.0);
        },
        std::invalid_argument
    );

    EXPECT_THROW(
        {
            SlaterPrimitive sp3(1, 0, 0, -1.0);
        },
        std::invalid_argument
    );
}

TEST(SlaterPrimitive, SettersGetters) {
    SlaterPrimitive sp(1, 0, 0, 1.0);

    sp.set_n(2);
    sp.set_l(1);
    sp.set_m(1);
    sp.set_alpha(3.0);

    EXPECT_EQ(sp.n(), 2);
    EXPECT_EQ(sp.l(), 1);
    EXPECT_EQ(sp.m(), 1);
    EXPECT_EQ(sp.alpha(), 3.0);
    EXPECT_NEAR(sp.normalization(), 18.0, 1e-6);

    // n = 2, l = 1, m = 1, alpha = 3.0
    // n < 1 -> throw
    EXPECT_THROW(
        {
            sp.set_n(0);
        },
        std::invalid_argument
    );

    // If n is set to 1, we have n = l = 1 which is invalid
    EXPECT_THROW(
        {
            sp.set_n(1);
        },
        std::invalid_argument
    );

    EXPECT_THROW(
        {
            sp.set_l(2);
        },
        std::invalid_argument
    );

    EXPECT_THROW(
        {
            sp.set_l(0);
        },
        std::invalid_argument
    );

    EXPECT_THROW(
        {
            sp.set_m(2);
        },
        std::invalid_argument
    );

    EXPECT_THROW(
        {
            sp.set_alpha(-1.0);
        },
        std::invalid_argument
    );
}