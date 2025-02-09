#include "BasisSet/contracted_orbital.hpp"
#include "BasisSet/contracted_orbital.interface.hpp"
#include <gtest/gtest.h>

TEST(ContractedSlater, Overlap) {}

TEST(ContractedSlater, Laplacian) {}

TEST(ContractedSlater, ElectronElectron) {
    ContractedSlater orbital1;
    orbital1.add_primitive(1, 1, 0, 0, 2);

    EXPECT_NEAR(
        electron_electron_integral(orbital1, orbital1, orbital1, orbital1),
        1.25, 1e-5);
}