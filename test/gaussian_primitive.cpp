#include "Orbitals/gaussian_primitive.hpp"
#include "./gaussian_tabulated_values.hpp"
#include "Eigen/src/Core/Matrix.h"
#include <gtest/gtest.h>

// @note We took the values from the following link:
// https://github.com/mortele/HartreeFock/blob/master/tests/integraltester.cpp

const double minimal_err = 1e-15;
const double relative_err = 1e-2;

TEST(OverlapIntegral, GaussianPrimitive) {
    auto cgp_1 = std::vector<GaussianPrimitive>();
    cgp_1.push_back(
        GaussianPrimitive(0, 0, 0, 1.00, Eigen::Vector3d{0.0, -1.5, 1.2}));
    cgp_1.push_back(
        GaussianPrimitive(1, 0, 0, 1.10, Eigen::Vector3d{0.1, -1.2, 2.5}));
    cgp_1.push_back(
        GaussianPrimitive(0, 1, 0, 0.30, Eigen::Vector3d{0.2, -1.0, -0.3}));
    cgp_1.push_back(
        GaussianPrimitive(0, 0, 1, 0.90, Eigen::Vector3d{0.3, -0.8, 0.1}));
    cgp_1.push_back(
        GaussianPrimitive(0, 0, 2, 2.20, Eigen::Vector3d{0.4, -0.6, -3.1}));
    cgp_1.push_back(
        GaussianPrimitive(0, 2, 0, 2.40, Eigen::Vector3d{0.5, -0.4, 3.8}));
    cgp_1.push_back(
        GaussianPrimitive(2, 0, 0, 3.10, Eigen::Vector3d{0.6, -0.2, 1.3}));
    cgp_1.push_back(
        GaussianPrimitive(1, 1, 0, 3.70, Eigen::Vector3d{0.7, 0.0, 2.4}));
    cgp_1.push_back(
        GaussianPrimitive(0, 1, 1, 1.30, Eigen::Vector3d{0.8, 0.2, 5.3}));
    cgp_1.push_back(
        GaussianPrimitive(1, 0, 1, 4.30, Eigen::Vector3d{0.9, 0.4, 1.2}));
    cgp_1.push_back(
        GaussianPrimitive(0, 0, 0, 0.03, Eigen::Vector3d{0.0, -1.5, 1.2}));

    auto m_exactOverlap = setupExactOverlapMatrix();

    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            double abs_err = std::max(
                std::abs(m_exactOverlap(i, j) * relative_err), minimal_err);
            EXPECT_NEAR(overlap_integral(cgp_1[i], cgp_1[j]),
                        m_exactOverlap(i, j), abs_err);
        }
    }
}

TEST(LaplacianIntegral, GaussianPrimitive) {
    auto cgp_1 = std::vector<GaussianPrimitive>();
    cgp_1.push_back(
        GaussianPrimitive(0, 0, 0, 1.00, Eigen::Vector3d{0.0, -1.5, 1.2}));
    cgp_1.push_back(
        GaussianPrimitive(1, 0, 0, 1.10, Eigen::Vector3d{0.1, -1.2, 2.5}));
    cgp_1.push_back(
        GaussianPrimitive(0, 1, 0, 0.30, Eigen::Vector3d{0.2, -1.0, -0.3}));
    cgp_1.push_back(
        GaussianPrimitive(0, 0, 1, 0.90, Eigen::Vector3d{0.3, -0.8, 0.1}));
    cgp_1.push_back(
        GaussianPrimitive(0, 0, 2, 2.20, Eigen::Vector3d{0.4, -0.6, -3.1}));
    cgp_1.push_back(
        GaussianPrimitive(0, 2, 0, 2.40, Eigen::Vector3d{0.5, -0.4, 3.8}));
    cgp_1.push_back(
        GaussianPrimitive(2, 0, 0, 3.10, Eigen::Vector3d{0.6, -0.2, 1.3}));
    cgp_1.push_back(
        GaussianPrimitive(1, 1, 0, 3.70, Eigen::Vector3d{0.7, 0.0, 2.4}));
    cgp_1.push_back(
        GaussianPrimitive(0, 1, 1, 1.30, Eigen::Vector3d{0.8, 0.2, 5.3}));
    cgp_1.push_back(
        GaussianPrimitive(1, 0, 1, 4.30, Eigen::Vector3d{0.9, 0.4, 1.2}));
    cgp_1.push_back(
        GaussianPrimitive(0, 0, 0, 0.03, Eigen::Vector3d{0.0, -1.5, 1.2}));

    auto m_exactKinetic = setupExactKineticMatrix();

    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            double abs_err = std::max(
                std::abs(m_exactKinetic(i, j) * relative_err), minimal_err);
            EXPECT_NEAR(-0.5 * laplacian_integral(cgp_1[i], cgp_1[j]),
                        m_exactKinetic(i, j), abs_err);
        }
    }
}

// TEST(ElectronNucleus, GaussianPrimitive) {
//     auto cgp_1 = std::vector<GaussianPrimitive>();
//     cgp_1.push_back(
//         GaussianPrimitive(0, 0, 0, 1.00, Eigen::Vector3d{0.0, -1.5, 1.2}));
//     cgp_1.push_back(
//         GaussianPrimitive(1, 0, 0, 1.10, Eigen::Vector3d{0.1, -1.2, 2.5}));
//     cgp_1.push_back(
//         GaussianPrimitive(0, 1, 0, 0.30, Eigen::Vector3d{0.2, -1.0, -0.3}));
//     cgp_1.push_back(
//         GaussianPrimitive(0, 0, 1, 0.90, Eigen::Vector3d{0.3, -0.8, 0.1}));
//     cgp_1.push_back(
//         GaussianPrimitive(0, 0, 2, 2.20, Eigen::Vector3d{0.4, -0.6, -3.1}));
//     cgp_1.push_back(
//         GaussianPrimitive(0, 2, 0, 2.40, Eigen::Vector3d{0.5, -0.4, 3.8}));
//     cgp_1.push_back(
//         GaussianPrimitive(2, 0, 0, 3.10, Eigen::Vector3d{0.6, -0.2, 1.3}));
//     cgp_1.push_back(
//         GaussianPrimitive(1, 1, 0, 3.70, Eigen::Vector3d{0.7, 0.0, 2.4}));
//     cgp_1.push_back(
//         GaussianPrimitive(0, 1, 1, 1.30, Eigen::Vector3d{0.8, 0.2, 5.3}));
//     cgp_1.push_back(
//         GaussianPrimitive(1, 0, 1, 4.30, Eigen::Vector3d{0.9, 0.4, 1.2}));
//     cgp_1.push_back(
//         GaussianPrimitive(0, 0, 0, 0.03, Eigen::Vector3d{0.0, -1.5, 1.2}));

//     auto m_exactElecNuc = setupExactElectronNucleusMatrix();

//     for (int i = 0; i < 10; i++) {
//         for (int j = 0; j < 10; j++) {
//             int atom =
//                 (359 * i + 295 * j - 42 * i * j + 120 * i * i -
//                  38 * i * j * i) %
//                 9; // some kind of weird black magic from the original code
//             atom = std::abs(atom);
//             auto position = cgp_1[atom].position();

//             double abs_err =
//                 std::max(std::abs(m_exactElecNuc(i, j) * 1e-1), minimal_err);
//             std::cout << "i: " << i << " j: " << j << " atom: " << atom
//                       << std::endl;
//             EXPECT_NEAR(electron_nucleus_integral(cgp_1[i], cgp_1[j],
//             position),
//                         m_exactElecNuc(i, j), abs_err);
//         }
//     }
// }

TEST(ElectronElectron, GaussianPrimitive) {
    auto m_elecElec = setupExactElectronElectronVector();

    for (int testIndex = 1; testIndex < 10; testIndex++) {
        double integral = 0;
        double difference = 0;

        Eigen::Vector3d posA;
        Eigen::Vector3d posC;
        Eigen::Vector3d posD;
        Eigen::Vector3d posB;

        GaussianPrimitive primitiveA;
        GaussianPrimitive primitiveB;
        GaussianPrimitive primitiveC;
        GaussianPrimitive primitiveD;

        switch (testIndex) {
        case 1:
            posA = {-0.5, 0, 0};
            posB = {-0.5, 0, 0};
            posC = {-0.5, 0, 0};
            posD = {-0.5, 0, 0};

            primitiveA = GaussianPrimitive(0, 0, 0, 13.0077, posA);
            primitiveB = GaussianPrimitive(0, 0, 0, 13.0077, posB);
            primitiveC = GaussianPrimitive(0, 0, 0, 13.0077, posC);
            primitiveD = GaussianPrimitive(0, 0, 0, 13.0077, posD);
            break;

        case 2:
            posA = {0.5, 0, 0};
            posB = {-0.5, 0, 0};
            posC = {-0.5, 0, 0};
            posD = {0.5, 0, 0};

            primitiveA = GaussianPrimitive(0, 0, 0, 13.0077, posA);
            primitiveB = GaussianPrimitive(0, 0, 0, 0.121949, posB);
            primitiveC = GaussianPrimitive(0, 0, 0, 0.444529, posC);
            primitiveD = GaussianPrimitive(0, 0, 0, 13.0077, posD);
            break;

        case 3:
            posA = {0.5, 0, 0};
            posB = {-0.5, 0, 0};
            posC = {-0.5, 0, 0};
            posD = {0.5, 0, 0};

            primitiveA = GaussianPrimitive(0, 0, 0, 13.0077, posA);
            primitiveB = GaussianPrimitive(0, 1, 0, 0.121949, posB);
            primitiveC = GaussianPrimitive(0, 1, 0, 0.444529, posC);
            primitiveD = GaussianPrimitive(0, 0, 0, 13.0077, posD);
            break;

        case 4:
            posA = {0.55, 1, 3};
            posB = {-0.52, 5, 6};
            posC = {-0.53, 1, 2};
            posD = {0.45, 2, 4};

            primitiveA = GaussianPrimitive(1, 0, 0, 13.0077, posA);
            primitiveB = GaussianPrimitive(0, 1, 0, 0.121949, posB);
            primitiveC = GaussianPrimitive(0, 1, 0, 0.444529, posC);
            primitiveD = GaussianPrimitive(0, 1, 0, 10.0077, posD);
            break;

        case 5:
            posA = {1.2, 2.3, 3.4};
            posB = {-1.3, 1.4, -2.4};
            posC = {2.3, 0.9, 3.2};
            posD = {5.0, 1.9, 1.2};

            primitiveA = GaussianPrimitive(0, 0, 0, 0.2, posA);
            primitiveB = GaussianPrimitive(0, 0, 0, 0.3, posB);
            primitiveC = GaussianPrimitive(0, 0, 0, 0.4, posC);
            primitiveD = GaussianPrimitive(0, 0, 0, 0.1, posD);
            break;

        case 6:
            posA = {1.2, 2.3, 3.4};
            posB = {-1.3, 1.4, -2.4};
            posC = {2.3, 0.9, 3.2};
            posD = {5.0, 1.9, 1.2};

            primitiveA = GaussianPrimitive(0, 0, 0, 0.2, posA);
            primitiveB = GaussianPrimitive(1, 0, 0, 0.3, posB);
            primitiveC = GaussianPrimitive(0, 0, 0, 0.4, posC);
            primitiveD = GaussianPrimitive(0, 0, 1, 0.1, posD);
            break;

        case 7:
            posA = {1.2, 2.3, 3.4};
            posB = {-1.3, 1.4, -2.4};
            posC = {2.3, 0.9, 3.2};
            posD = {5.0, 1.9, 1.2};

            primitiveA = GaussianPrimitive(0, 0, 0, 0.2, posA);
            primitiveB = GaussianPrimitive(1, 0, 0, 0.3, posB);
            primitiveC = GaussianPrimitive(0, 2, 0, 0.4, posC);
            primitiveD = GaussianPrimitive(0, 0, 1, 0.1, posD);
            break;

        case 8:
            posA = {1.2, 2.3, 3.4};
            posB = {-1.3, 1.4, -2.4};
            posC = {2.3, 0.9, 3.2};
            posD = {5.0, 1.9, 1.2};

            primitiveA = GaussianPrimitive(1, 1, 0, 0.2, posA);
            primitiveB = GaussianPrimitive(2, 0, 0, 0.3, posB);
            primitiveC = GaussianPrimitive(2, 0, 0, 0.4, posC);
            primitiveD = GaussianPrimitive(2, 0, 0, 0.1, posD);
            break;

        case 9:
            posA = {1.2, 2.3, 3.4};
            posB = {-1.3, 1.4, -2.4};
            posC = {2.3, 0.9, 3.2};
            posD = {5.0, 1.9, 1.2};

            primitiveA = GaussianPrimitive(1, 1, 0, 0.2, posA);
            primitiveB = GaussianPrimitive(0, 2, 0, 0.3, posB);
            primitiveC = GaussianPrimitive(0, 2, 0, 0.4, posC);
            primitiveD = GaussianPrimitive(2, 0, 0, 0.1, posD);
            break;

        default:
            break;
        }

        double abs_err = std::max(
            std::abs(electron_electron_integral(primitiveA, primitiveB,
                                                primitiveC, primitiveD) *
                     relative_err),
            minimal_err);
        EXPECT_NEAR(electron_electron_integral(primitiveA, primitiveB,
                                               primitiveC, primitiveD),
                    m_elecElec(testIndex), abs_err);
    }
}