#include <gtest/gtest.h>
#include "Integrators/laplacian_integral.hpp"

TEST(LaplacianIntegral, SlaterPrimitiveZero) {
    SlaterPrimitive orbital1(1, 0, 0, 0.5);
    SlaterPrimitive orbital2(1, 1, 0, 0.5);

    EXPECT_NEAR(laplacian_integral(orbital1, orbital2), 0., 1e-6);

    orbital1.set_l(1); orbital1.set_m(1);
    EXPECT_NEAR(laplacian_integral(orbital1, orbital2), 0., 1e-6);
}

TEST(LaplacianIntegral, SlaterPrimitiveValue) {
    SlaterPrimitive orbital1(1, 0, 0, 0.5);
    SlaterPrimitive orbital2(5, 0, 0, 1.3);

    EXPECT_NEAR(laplacian_integral(orbital1, orbital2), -0.588023882237, 1e-6);

    orbital1.set_n(5); orbital2.set_n(3);
    orbital1.set_l(1); orbital2.set_l(1);
    EXPECT_NEAR(laplacian_integral(orbital1, orbital2), 33.8295221138913, 1e-6);

    orbital1.set_m(1); orbital2.set_m(1);
    EXPECT_NEAR(laplacian_integral(orbital1, orbital2), 33.8295221138913, 1e-6);

    orbital1.set_l(2); orbital2.set_l(2);
    EXPECT_NEAR(laplacian_integral(orbital1, orbital2), -13.212388465, 1e-6);
}