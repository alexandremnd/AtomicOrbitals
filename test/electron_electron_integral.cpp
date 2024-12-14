#include <gtest/gtest.h>
#include "BasisSet/slater_primitive.hpp"
#include "Integrators/electron_electron_integral.hpp"

TEST(ElectronElectronIntegralTest, SlaterPrimitive) {
    SlaterPrimitive orbital1(3, 2, 0, 1.5);
    SlaterPrimitive orbital2(3, 2, 0, 1);
    SlaterPrimitive orbital3(3, 2, 0, 0.3);
    SlaterPrimitive orbital4(3, 2, 0, 1);

    double result = electron_electron_integral(orbital1, orbital2, orbital3, orbital4);
    EXPECT_NEAR(result, 17.2833787319263, 1e-5);
}