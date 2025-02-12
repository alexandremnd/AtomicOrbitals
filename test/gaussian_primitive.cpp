#include "BasisSet/gaussian_primitive.hpp"
#include "BasisSet/slater_primitive.hpp"
#include "Eigen/src/Core/Matrix.h"
#include <gtest/gtest.h>

TEST(OverlapIntegralTest, GaussianPrimitive) {
    auto gp_1 = GaussianPrimitive(0, 0, 1, 1.0, Eigen::Vector3d(0, 0, 0));
    auto gp_2 = GaussianPrimitive(0, 0, 1, 2.0, Eigen::Vector3d(1, 1, 1));

    EXPECT_NEAR(overlap_integral(gp_1, gp_2) / gp_1.normalization() /
                    gp_2.normalization(),
                -0.00805715016827179, 1e-6);

    auto gp_3 = GaussianPrimitive(1, 2, 4, 0.5, Eigen::Vector3d(-1, -5, -4));
    auto gp_4 = GaussianPrimitive(3, 3, 5, 0.5, Eigen::Vector3d(3, 3, 3));
    EXPECT_NEAR(overlap_integral(gp_3, gp_4) / gp_3.normalization() /
                    gp_4.normalization(),
                -0.0000548547, 1e-6);
}

TEST(LaplacianIntegral, GaussianPrimitive) {
    auto gp_1 = GaussianPrimitive(0, 0, 1, 1.0, Eigen::Vector3d(0, 0, 0));
    auto gp_2 = GaussianPrimitive(0, 0, 1, 2.0, Eigen::Vector3d(1, 1, 1));

    EXPECT_NEAR(laplacian_integral(gp_1, gp_2) / gp_1.normalization() /
                    gp_2.normalization(),
                0.0966858, 1e-6);

    auto gp_4 = GaussianPrimitive(2, 0, 1, 1.0, Eigen::Vector3d(0, 0, 0));
    auto gp_3 = GaussianPrimitive(0, 2, 1, 2.0, Eigen::Vector3d(1, -1, 1));

    EXPECT_NEAR(laplacian_integral(gp_4, gp_3) / gp_4.normalization() /
                    gp_3.normalization(),
                0.021552, 1e-6);

    auto gp_5 = GaussianPrimitive(2, 3, 1, 0.7, Eigen::Vector3d(1, 0, 0));
    auto gp_6 = GaussianPrimitive(0, 3, 1, 0.5, Eigen::Vector3d(1, -1, 1));

    EXPECT_NEAR(laplacian_integral(gp_5, gp_6) / gp_5.normalization() /
                    gp_6.normalization(),
                0.076246, 1e-6);
}

TEST(ElectronNucleus, GaussianPrimitive) {
    auto gp_1 = GaussianPrimitive(0, 0, 1, 1.0, Eigen::Vector3d(0, 0, 0));
    auto gp_2 = GaussianPrimitive(0, 0, 1, 2.0, Eigen::Vector3d(0, 0, 0));

    EXPECT_NEAR(electron_nucleus_integral(gp_1, gp_2, Eigen::Vector3d::Zero()) /
                    gp_1.normalization() / gp_2.normalization(),
                0.00473852, 1e-6);

    auto gp_4 = GaussianPrimitive(2, 0, 1, 1.0, Eigen::Vector3d(0, 0, 0));
    auto gp_3 = GaussianPrimitive(0, 2, 1, 2.0, Eigen::Vector3d(1, -1, 1));

    EXPECT_NEAR(laplacian_integral(gp_4, gp_3) / gp_4.normalization() /
                    gp_3.normalization(),
                0.021552, 1e-6);

    auto gp_5 = GaussianPrimitive(2, 3, 1, 0.7, Eigen::Vector3d(1, 0, 0));
    auto gp_6 = GaussianPrimitive(0, 3, 1, 0.5, Eigen::Vector3d(1, -1, 1));

    EXPECT_NEAR(laplacian_integral(gp_5, gp_6) / gp_5.normalization() /
                    gp_6.normalization(),
                0.076246, 1e-6);
}