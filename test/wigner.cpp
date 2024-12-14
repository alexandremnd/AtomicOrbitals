#include <gtest/gtest.h>
#include "Maths/wigner_3j.hpp"

TEST(Wigner3jTest, HandlesZeroCase) {
    double r1 = Math::wigner_3j(1, 1, 1, 1, 0, 0);
    double r2 = Math::wigner_3j(1, 1, 1, 0, 1, 0);
    double r3 = Math::wigner_3j(1, 1, 1, 0, 0, 1);

    double r4 = Math::wigner_3j(1, 1, 1, 2, -1, -1);
    double r5 = Math::wigner_3j(1, 1, 1, -1, 2, -1);
    double r6 = Math::wigner_3j(1, 1, 1, -1, -1, 2);

    double r7 = Math::wigner_3j(3, 1, 1, 0, 0, 0);
    double r8 = Math::wigner_3j(1, 3, 1, 0, 0, 0);
    double r9 = Math::wigner_3j(1, 1, 3, 0, 0, 0);

    ASSERT_NEAR(r1, 0, 1e-6);
    ASSERT_NEAR(r2, 0, 1e-6);
    ASSERT_NEAR(r3, 0, 1e-6);

    ASSERT_NEAR(r4, 0, 1e-6);
    ASSERT_NEAR(r5, 0, 1e-6);
    ASSERT_NEAR(r6, 0, 1e-6);

    ASSERT_NEAR(r7, 0, 1e-6);
    ASSERT_NEAR(r8, 0, 1e-6);
    ASSERT_NEAR(r9, 0, 1e-6);
}

TEST(Wigner3jTest, CheckReturnedValue) {
    double r1 = Math::wigner_3j(1, 1, 1, 1, -1, 0);
    double r2 = Math::wigner_3j(2, 3, 4, 1, -2, 1);
    double r3 = Math::wigner_3j(2, 3, 3, -1, 3, -2);

    ASSERT_NEAR(r1, 0.408248, 1e-6);
    ASSERT_NEAR(r2, 0.197203, 1e-6);
    ASSERT_NEAR(r3, -0.243975, 1e-6);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}