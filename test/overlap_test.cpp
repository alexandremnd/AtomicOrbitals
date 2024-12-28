#include "gtest/gtest.h"
#include <gtest/gtest.h>
#include "BasisSet/slater_primitive.hpp"
#include "Integrators/overlap_integral.hpp"

TEST(OverlapIntegralTest, SlaterPrimitiveZeroCase) {
    SlaterPrimitive orbital1(2, 0, 0, 0.5);
    SlaterPrimitive orbital2(5, 1, 0, 1.3);

    EXPECT_NEAR(overlap_integral(orbital1, orbital2), 0., 1e-6);

    orbital1.set_l(1); orbital1.set_m(0);
    orbital2.set_l(1); orbital2.set_m(1);

    EXPECT_NEAR(overlap_integral(orbital1, orbital2), 0., 1e-6);
}

TEST(OverlapIntegralTest, SlaterPrimitiveValue) {
    SlaterPrimitive orbital1(1, 0, 0, 0.5);
    SlaterPrimitive orbital2(5, 0, 0, 1.3);

    EXPECT_NEAR(overlap_integral(orbital1, orbital2), 11.76047764474325, 1e-5);

    orbital1.set_n(3); orbital2.set_n(8);
    EXPECT_NEAR(overlap_integral(orbital1, orbital2), 34505.28801422153, 1e-3);
}

TEST(OverlapIntegralTest, SlaterPrimitiveOffsetThrow) {
    SlaterPrimitive orbital1(1, 0, 0, 0.5);
    SlaterPrimitive orbital2(1, 0, 0, 1.3);

    EXPECT_DEATH({
        overlap_integral(orbital1, orbital2, -3);
    }, ".*");

    EXPECT_EXIT({
        overlap_integral(orbital1, orbital2, -2);
        exit(0);
    }, ::testing::ExitedWithCode(0), ".*");

    EXPECT_EXIT({
        overlap_integral(orbital1, orbital2, -1);
        exit(0);
    }, ::testing::ExitedWithCode(0), ".*");

    EXPECT_EXIT({
        overlap_integral(orbital1, orbital2, 0);
        exit(0);
    }, ::testing::ExitedWithCode(0), ".*");


    SlaterPrimitive orbital3(2, 0, 0, 0.5);
    SlaterPrimitive orbital4(1, 0, 0, 0.5);

    EXPECT_DEATH(overlap_integral(orbital3, orbital4, -4), ".*");

    EXPECT_EXIT({
        overlap_integral(orbital3, orbital4, -3);
        exit(0);
    }, ::testing::ExitedWithCode(0), ".*");

    EXPECT_EXIT({
        overlap_integral(orbital3, orbital4, -2);
        exit(0);
    }, ::testing::ExitedWithCode(0), ".*");

    EXPECT_EXIT({
        overlap_integral(orbital3, orbital4, -1);
        exit(0);
    }, ::testing::ExitedWithCode(0), ".*");

    EXPECT_EXIT({
        overlap_integral(orbital3, orbital4, 0);
        exit(0);
    }, ::testing::ExitedWithCode(0), ".*");
}