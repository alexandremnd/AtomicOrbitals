#include <gtest/gtest.h>
#include "BasisSet/slater_primitive.hpp"
#include "Integrators/slater_braket.hpp"

// For floating point comparison, we allow a relative error of around 1e-6
TEST(SlaterBraketTest, RadialIntegral) {
    SlaterPrimitive orbital1(1, 1, 1, 1);
    SlaterPrimitive orbital2(3, 1, 1, 1);
    SlaterPrimitive orbital3(2, 1, 1, 1);
    SlaterPrimitive orbital4(2, 1, 1, 1);

    double result = radial_integral(orbital1, orbital2, orbital3, orbital4, 0);
    EXPECT_NEAR(result, 0.246826, 1e-6);

    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 2);
    EXPECT_NEAR(result, 0.106201, 1e-6);

    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 3);
    EXPECT_NEAR(result, 0.0807018, 1e-6);

    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 4);
    EXPECT_NEAR(result, 0.0646995, 1e-6);


    orbital1.set_n(2);
    orbital3.set_n(1);
    orbital2.set_n(1);
    orbital4.set_n(4);
    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 0);
    EXPECT_NEAR(result, 0.246826, 1e-6);

    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 2);
    EXPECT_NEAR(result, 0.106201, 1e-6);

    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 3);
    EXPECT_NEAR(result, 0.0807018, 1e-6);

    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 4);
    EXPECT_NEAR(result, 0.0646995, 1e-6);

    orbital1 = SlaterPrimitive(1, 1, 1, 0.2);
    orbital3 = SlaterPrimitive(4, 1, 1, 0.5);

    orbital2 = SlaterPrimitive(6, 1, 1, 0.1);
    orbital4 = SlaterPrimitive(2, 1, 1, 0.8);

    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 0);
    EXPECT_NEAR(result, 1.022197954433982e7, 5);

    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 2);
    EXPECT_NEAR(result, 5.289151185199601e6, 5);

    result = radial_integral(orbital1, orbital2, orbital3, orbital4, 5);
    EXPECT_NEAR(result, 2.894457084101954e6, 5);
}

TEST(SlaterBraketTest, AngularIntegralTest) {
    SlaterPrimitive orbital1(1, 0, 0, 1.5);
    SlaterPrimitive orbital2(3, 0, 0, 0.3);

    double result = angular_integral(orbital1, orbital2, 0, 0);
    EXPECT_NEAR(result, 0.282095, 1e-6);



    orbital1 = SlaterPrimitive(1, 2, 0, 1.5);
    orbital2 = SlaterPrimitive(3, 0, 0, 0.3);

    result = angular_integral(orbital1, orbital2, 0, 0);
    EXPECT_NEAR(result, 0.0, 1e-6);

    result = angular_integral(orbital1, orbital2, 2, 0);
    EXPECT_NEAR(result, 0.282095, 1e-6);



    orbital1 = SlaterPrimitive(1, 2, 0, 1.5);
    orbital2 = SlaterPrimitive(3, 1, 0, 0.3);

    result = angular_integral(orbital1, orbital2, 0, 0);
    EXPECT_NEAR(result, 0.0, 1e-6);

    result = angular_integral(orbital1, orbital2, 1, 0);
    EXPECT_NEAR(result, 0.252313, 1e-6);



    orbital1 = SlaterPrimitive(1, 3, -2, 1.5);
    orbital2 = SlaterPrimitive(3,4, 2, 0.3);

    result = angular_integral(orbital1, orbital2, 7, -4);
    EXPECT_NEAR(result, 0.25085, 1e-6);

    result = angular_integral(orbital1, orbital2, 5, 0);
    EXPECT_NEAR(result, 0.0, 1e-6);



    orbital1 = SlaterPrimitive(1, 3, -2, 1.5);
    orbital2 = SlaterPrimitive(3,4, 3, 0.3);

    result = angular_integral(orbital1, orbital2, 5, -5);
    EXPECT_NEAR(result, -0.212007, 1e-6);

    orbital1 = SlaterPrimitive(1, 3, -2, 1.5);
    orbital2 = SlaterPrimitive(3,3, 2, 0.3);

    result = angular_integral(orbital1, orbital2, 4, -4);
    EXPECT_NEAR(result, 0.214561, 1e-6);

}