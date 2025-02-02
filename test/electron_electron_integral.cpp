#include <gtest/gtest.h>
#include "BasisSet/slater_primitive.hpp"
#include "Integrators/electron_electron_integral.hpp"

TEST(ElectronElectronIntegralTest, SlaterPrimitive) {
    SlaterPrimitive orbital1(3, 2, 0, 1.5);
    SlaterPrimitive orbital2(3, 2, 0, 1);
    SlaterPrimitive orbital3(3, 2, 0, 0.3);
    SlaterPrimitive orbital4(3, 2, 0, 1);

    double normalization = orbital1.normalization() * orbital2.normalization() * orbital3.normalization() * orbital4.normalization();
    double result = electron_electron_integral(orbital1, orbital2, orbital3, orbital4);
    EXPECT_NEAR(result / normalization, 17.2833787319263, 1e-5);


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
    normalization = orbital1.normalization() * orbital2.normalization() * orbital3.normalization() * orbital4.normalization();
    result = electron_electron_integral(orbital1, orbital2, orbital3, orbital4);
    EXPECT_NEAR(result / normalization,  0.00043001174698042, 1e-7);

    orbital1 = SlaterPrimitive(1, 0, 0, 2);
    orbital2 = SlaterPrimitive(1, 0, 0, 2);
    orbital3 = SlaterPrimitive(1, 0, 0, 2);
    orbital4 = SlaterPrimitive(1, 0, 0, 2);
    result = electron_electron_integral(orbital1, orbital2, orbital3, orbital4);
    EXPECT_NEAR(result, 1.25, 1e-6);


    orbital1 = SlaterPrimitive(2, 1, 0, 2);
    orbital2 = SlaterPrimitive(3, 2, 0, 2);
    orbital3 = SlaterPrimitive(4, 2, -1, 2);
    orbital4 = SlaterPrimitive(2, 1, -1, 2);
    result = electron_electron_integral(orbital1, orbital2, orbital3, orbital4);
    EXPECT_NEAR(result, 0.0, 1e-7);
}