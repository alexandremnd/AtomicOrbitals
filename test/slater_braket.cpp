#include <gtest/gtest.h>
#include "BasisSet/slater_primitive.hpp"
#include "Integrators/slater_braket.hpp"

// For floating point comparison, we allow a relative error of around 1e-6
TEST(SlaterBraketTest, RadialIntegral) {
    SlaterPrimitive orbital1(2, 1, 1, 1);
    SlaterPrimitive orbital2(3, 1, 1, 1);
    SlaterPrimitive orbital3(2, 1, 1, 1);
    SlaterPrimitive orbital4(2, 1, 1, 1);

    double result = radial_integral(orbital1, orbital2, orbital3, orbital4, 0);
    EXPECT_NEAR(result, 0.459777832031250, 1e-6);

    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 1);
    EXPECT_NEAR(result, 0.304870605468750, 1e-6);

    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 2);
    EXPECT_NEAR(result, 0.222473144531250, 1e-6);


    orbital1 = SlaterPrimitive(2, 1, 1, 3);
    orbital2 = SlaterPrimitive(3, 1, 1, 6);
    orbital3 = SlaterPrimitive(3, 1, 1, 9);
    orbital4 = SlaterPrimitive(2, 1, 1, 12);

    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 0);
    EXPECT_NEAR(result, 3.00294126932456e-10, 1e-15);

    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 1);
    EXPECT_NEAR(result, 1.92980029768121e-10, 1e-15);

    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 2);
    EXPECT_NEAR(result, 1.37973923728128e-10, 1e-15);

    orbital1 = SlaterPrimitive(4, 3, 0, 3);
    orbital2 = SlaterPrimitive(3, 2, 0, 6);
    orbital3 = SlaterPrimitive(5, 4, 0, 9);
    orbital4 = SlaterPrimitive(3, 2, 0, 12);

    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 0);
    EXPECT_NEAR(result, 9.01607666199089e-12, 1e-18);

    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 1);
    EXPECT_NEAR(result, 4.79826773710592e-12, 1e-18);

    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 2);
    EXPECT_NEAR(result, 2.98402484593426e-12, 1e-18);

    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 3);
    EXPECT_NEAR(result, 2.07510878703165e-12, 1e-18);

    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 4);
    EXPECT_NEAR(result, 1.56099109151659e-12, 1e-18);


    orbital1 = SlaterPrimitive(2, 1, 0, 50);
    orbital2 = SlaterPrimitive(3, 1, 0, 47);
    orbital3 = SlaterPrimitive(2, 1, 0, 50);
    orbital4 = SlaterPrimitive(3, 1, 0, 47);

    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 0);
    EXPECT_NEAR(result, 3.74105950679226e-19, 1e-25);

    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 1);
    EXPECT_NEAR(result, 2.36867167366852e-19, 1e-25);

    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 2);
    EXPECT_NEAR(result, 1.68203740085704e-19, 1e-25);
}

TEST(SlaterBraketTest, AngularIntegralTest) {
    SlaterPrimitive orbital1(1, 0, 0, 1.5);
    SlaterPrimitive orbital2(3, 0, 0, 0.3);

    double result = angular_integral(orbital1, orbital2, 0, 0);
    EXPECT_NEAR(result, 0.282095, 1e-6);



    orbital1 = SlaterPrimitive(3, 2, 0, 1.5);
    orbital2 = SlaterPrimitive(3, 0, 0, 0.3);

    result = angular_integral(orbital1, orbital2, 0, 0);
    EXPECT_NEAR(result, 0.0, 1e-6);

    result = angular_integral(orbital1, orbital2, 2, 0);
    EXPECT_NEAR(result, 0.282094791773878, 1e-6);



    orbital1 = SlaterPrimitive(3, 2, 0, 1.5);
    orbital2 = SlaterPrimitive(3, 1, 0, 0.3);

    result = angular_integral(orbital1, orbital2, 0, 0);
    EXPECT_NEAR(result, 0.0, 1e-6);

    result = angular_integral(orbital1, orbital2, 1, 0);
    EXPECT_NEAR(result, 0.252313, 1e-6);



    orbital1 = SlaterPrimitive(4, 3, -2, 1.5);
    orbital2 = SlaterPrimitive(5,4, 2, 0.3);

    result = angular_integral(orbital1, orbital2, 7, -4);
    EXPECT_NEAR(result, 0.25085, 1e-6);

    result = angular_integral(orbital1, orbital2, 5, 0);
    EXPECT_NEAR(result, 0.0, 1e-6);



    orbital1 = SlaterPrimitive(4, 3, -2, 1.5);
    orbital2 = SlaterPrimitive(5,4, 3, 0.3);

    result = angular_integral(orbital1, orbital2, 5, -5);
    EXPECT_NEAR(result, -0.212007, 1e-6);

    orbital1 = SlaterPrimitive(4, 3, -2, 1.5);
    orbital2 = SlaterPrimitive(4,3, 2, 0.3);

    result = angular_integral(orbital1, orbital2, 4, -4);
    EXPECT_NEAR(result, 0.214561, 1e-6);

}

TEST(SlaterBraketTest, Radial_LNegative) {
    SlaterPrimitive orbital1(2, 1, 1, 1);
    SlaterPrimitive orbital2(3, 1, 1, 1);
    SlaterPrimitive orbital3(2, 1, 1, 1);
    SlaterPrimitive orbital4(2, 1, 1, 1);

    EXPECT_DEATH(radial_integral(orbital1, orbital2, orbital3, orbital4, -1), ".*");

    EXPECT_EXIT({
        radial_integral(orbital1, orbital2, orbital3, orbital4, 0);
        exit(0);
    }, ::testing::ExitedWithCode(0), ".*");

    EXPECT_EXIT({
        radial_integral(orbital1, orbital2, orbital3, orbital4, 1);
        exit(0);
    }, ::testing::ExitedWithCode(0), ".*");
}

TEST(SlaterBraketTest, Angular_LNegative) {
    SlaterPrimitive orbital1(2, 1, 1, 1);
    SlaterPrimitive orbital2(3, 1, 1, 1);
    SlaterPrimitive orbital3(2, 1, 1, 1);
    SlaterPrimitive orbital4(2, 1, 1, 1);

    EXPECT_DEATH(angular_integral(orbital1, orbital2, -1, 0), ".*");

    EXPECT_EXIT({
        angular_integral(orbital1, orbital2, 0, 0);
        exit(0);
    }, ::testing::ExitedWithCode(0), ".*");

    EXPECT_EXIT({
        angular_integral(orbital1, orbital2, 1, 0);
        exit(0);
    }, ::testing::ExitedWithCode(0), ".*");
}

TEST(SlaterBraketTest, Radial_LTooLarge) {
    SlaterPrimitive orbital1(2, 1, 1, 1);
    SlaterPrimitive orbital2(3, 1, 1, 1);
    SlaterPrimitive orbital3(2, 1, 1, 1);
    SlaterPrimitive orbital4(2, 1, 1, 1);

    EXPECT_DEATH(radial_integral(orbital1, orbital2, orbital3, orbital4, 4), ".*");
    EXPECT_DEATH(radial_integral(orbital1, orbital2, orbital3, orbital4, 3), ".*");

    EXPECT_EXIT({
        radial_integral(orbital1, orbital2, orbital3, orbital4, 2);
        exit(0);
    }, ::testing::ExitedWithCode(0), ".*");

    EXPECT_EXIT({
        radial_integral(orbital1, orbital2, orbital3, orbital4, 1);
        exit(0);
    }, ::testing::ExitedWithCode(0), ".*");

    EXPECT_EXIT({
        radial_integral(orbital1, orbital2, orbital3, orbital4, 0);
        exit(0);
    }, ::testing::ExitedWithCode(0), ".*");
}


TEST(SlaterBraketTest, Angular_MTooLarge) {
    SlaterPrimitive orbital1(2, 1, 1, 1);
    SlaterPrimitive orbital2(3, 1, 1, 1);

    EXPECT_DEATH(angular_integral(orbital1, orbital2, 1, 2), ".*");

    EXPECT_EXIT({
        angular_integral(orbital1, orbital2, 0, 0);
        exit(0);
    }, ::testing::ExitedWithCode(0), ".*");

    EXPECT_EXIT({
        angular_integral(orbital1, orbital2, 1, 0);
        exit(0);
    }, ::testing::ExitedWithCode(0), ".*");

    EXPECT_EXIT({
        angular_integral(orbital1, orbital2, 1, -1);
        exit(0);
    }, ::testing::ExitedWithCode(0), ".*");

    EXPECT_EXIT({
        angular_integral(orbital1, orbital2, 1, 1);
        exit(0);
    }, ::testing::ExitedWithCode(0), ".*");
}