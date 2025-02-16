#include <BasisSet/slater_primitive.hpp>
#include <gtest/gtest.h>
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

    EXPECT_THROW({ SlaterPrimitive sp3(0, 0, 0, 1.0); }, std::invalid_argument);

    EXPECT_THROW({ SlaterPrimitive sp3(1, 2, 0, 1.0); }, std::invalid_argument);

    EXPECT_THROW({ SlaterPrimitive sp3(1, 0, 3, 1.0); }, std::invalid_argument);

    EXPECT_THROW(
        { SlaterPrimitive sp3(1, 0, 0, -1.0); }, std::invalid_argument);
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
    EXPECT_THROW({ sp.set_n(0); }, std::invalid_argument);

    // If n is set to 1, we have n = l = 1 which is invalid
    EXPECT_THROW({ sp.set_n(1); }, std::invalid_argument);

    EXPECT_THROW({ sp.set_l(2); }, std::invalid_argument);

    EXPECT_THROW({ sp.set_l(0); }, std::invalid_argument);

    EXPECT_THROW({ sp.set_m(2); }, std::invalid_argument);

    EXPECT_THROW({ sp.set_alpha(-1.0); }, std::invalid_argument);
}

TEST(OverlapIntegralTest, SlaterPrimitiveZeroCase) {
    SlaterPrimitive orbital1(2, 0, 0, 0.5);
    SlaterPrimitive orbital2(5, 1, 0, 1.3);

    EXPECT_NEAR(overlap_integral(orbital1, orbital2), 0., 1e-6);

    orbital1.set_l(1);
    orbital1.set_m(0);
    orbital2.set_l(1);
    orbital2.set_m(1);

    EXPECT_NEAR(overlap_integral(orbital1, orbital2), 0., 1e-6);
}

TEST(OverlapIntegralTest, SlaterPrimitiveValue) {
    SlaterPrimitive orbital1(1, 0, 0, 0.5);
    SlaterPrimitive orbital2(5, 0, 0, 1.3);

    EXPECT_NEAR(overlap_integral(orbital1, orbital2), 11.76047764474325, 1e-5);

    orbital1.set_n(3);
    orbital2.set_n(8);
    EXPECT_NEAR(overlap_integral(orbital1, orbital2), 34505.28801422153, 1e-3);
}

TEST(OverlapIntegralTest, SlaterPrimitiveOffsetThrow) {
    SlaterPrimitive orbital1(1, 0, 0, 0.5);
    SlaterPrimitive orbital2(1, 0, 0, 1.3);

    EXPECT_DEATH({ overlap_integral(orbital1, orbital2, -3); }, ".*");

    EXPECT_EXIT(
        {
            overlap_integral(orbital1, orbital2, -2);
            exit(0);
        },
        ::testing::ExitedWithCode(0), ".*");

    EXPECT_EXIT(
        {
            overlap_integral(orbital1, orbital2, -1);
            exit(0);
        },
        ::testing::ExitedWithCode(0), ".*");

    EXPECT_EXIT(
        {
            overlap_integral(orbital1, orbital2, 0);
            exit(0);
        },
        ::testing::ExitedWithCode(0), ".*");

    SlaterPrimitive orbital3(2, 0, 0, 0.5);
    SlaterPrimitive orbital4(1, 0, 0, 0.5);

    EXPECT_DEATH(overlap_integral(orbital3, orbital4, -4), ".*");

    EXPECT_EXIT(
        {
            overlap_integral(orbital3, orbital4, -3);
            exit(0);
        },
        ::testing::ExitedWithCode(0), ".*");

    EXPECT_EXIT(
        {
            overlap_integral(orbital3, orbital4, -2);
            exit(0);
        },
        ::testing::ExitedWithCode(0), ".*");

    EXPECT_EXIT(
        {
            overlap_integral(orbital3, orbital4, -1);
            exit(0);
        },
        ::testing::ExitedWithCode(0), ".*");

    EXPECT_EXIT(
        {
            overlap_integral(orbital3, orbital4, 0);
            exit(0);
        },
        ::testing::ExitedWithCode(0), ".*");
}

TEST(LaplacianIntegral, SlaterPrimitiveZero) {
    SlaterPrimitive orbital1(2, 0, 0, 0.5);
    SlaterPrimitive orbital2(2, 1, 0, 0.5);

    EXPECT_NEAR(laplacian_integral(orbital1, orbital2), 0., 1e-6);

    orbital1.set_l(1);
    orbital1.set_m(1);
    EXPECT_NEAR(laplacian_integral(orbital1, orbital2), 0., 1e-6);
}

TEST(LaplacianIntegral, SlaterPrimitiveValue) {
    SlaterPrimitive orbital1(1, 0, 0, 0.5);
    SlaterPrimitive orbital2(5, 0, 0, 1.3);

    EXPECT_NEAR(laplacian_integral(orbital1, orbital2), -0.588023882237, 1e-6);

    orbital1.set_n(5);
    orbital2.set_n(3);
    orbital1.set_l(1);
    orbital2.set_l(1);
    EXPECT_NEAR(laplacian_integral(orbital1, orbital2), 33.8295221138913, 1e-6);

    orbital1.set_m(1);
    orbital2.set_m(1);
    EXPECT_NEAR(laplacian_integral(orbital1, orbital2), 33.8295221138913, 1e-6);

    orbital1.set_l(2);
    orbital2.set_l(2);
    EXPECT_NEAR(laplacian_integral(orbital1, orbital2), -13.212388465, 1e-6);
}

TEST(ElectronElectronIntegralTest, SlaterPrimitive) {
    SlaterPrimitive orbital1(3, 2, 0, 1.5);
    SlaterPrimitive orbital2(3, 2, 0, 1);
    SlaterPrimitive orbital3(3, 2, 0, 0.3);
    SlaterPrimitive orbital4(3, 2, 0, 1);
    double result =
        electron_electron_integral(orbital1, orbital2, orbital3, orbital4);
    EXPECT_NEAR(result, 17.2833787319263, 1e-5);

    orbital1 = SlaterPrimitive(2, 1, 0, 2);
    orbital2 = SlaterPrimitive(2, 1, 0, 2);
    orbital3 = SlaterPrimitive(2, 1, -1, 2);
    orbital4 = SlaterPrimitive(2, 1, -1, 2);
    result = electron_electron_integral(orbital1, orbital2, orbital3, orbital4);
    EXPECT_NEAR(result, 0.0, 1e-5);

    orbital1 = SlaterPrimitive(2, 1, 0, 2);
    orbital2 = SlaterPrimitive(2, 1, 0, 2);
    orbital3 = SlaterPrimitive(2, 1, 0, 2);
    orbital4 = SlaterPrimitive(2, 1, 0, 2);
    result = electron_electron_integral(orbital1, orbital2, orbital3, orbital4);
    EXPECT_NEAR(result, 0.00043001174698042, 1e-7);

    orbital1 = SlaterPrimitive(1, 0, 0, 2);
    orbital2 = SlaterPrimitive(1, 0, 0, 2);
    orbital3 = SlaterPrimitive(1, 0, 0, 2);
    orbital4 = SlaterPrimitive(1, 0, 0, 2);
    result = electron_electron_integral(orbital1, orbital2, orbital3, orbital4);
    double normalization = orbital1.normalization() * orbital2.normalization() *
                           orbital3.normalization() * orbital4.normalization();
    EXPECT_NEAR(result, 1.25 / normalization, 1e-6);

    orbital1 = SlaterPrimitive(2, 1, 0, 2);
    orbital2 = SlaterPrimitive(3, 2, 0, 2);
    orbital3 = SlaterPrimitive(4, 2, -1, 2);
    orbital4 = SlaterPrimitive(2, 1, -1, 2);
    result = electron_electron_integral(orbital1, orbital2, orbital3, orbital4);
    EXPECT_NEAR(result, 0.0, 1e-7);
}