#pragma once

#include "BasisSet/slater_primitive.hpp"
#include "BasisSet/gaussian_primitive.hpp"
#include "BasisSet/contracted_orbital.hpp"
#include "Eigen/Dense"

double electron_nucleus_integral(const GaussianPrimitive& orbital1,
                                const GaussianPrimitive& orbital2,
                                const Eigen::Vector3d& nucleus_position);

double electron_nucleus_integral(const SlaterPrimitive& orbital1,
                                const SlaterPrimitive& orbital2,
                                const Eigen::Vector3d& nucleus_position);

template <DerivedFromOrbital PrimitiveType>
double electron_nucleus_integral(const ContractedOrbital<PrimitiveType>& orbital1,
                                const ContractedOrbital<PrimitiveType>& orbital2,
                                const Eigen::Vector3d& nucleus_position) {
    double integral = 0.0;
    for (size_t i = 0; i < orbital1.size(); i++) {
        for (size_t j = 0; j < orbital2.size(); j++) {
            integral += orbital1.get_coefficient(i) * orbital2.get_coefficient(j) *
                        electron_nucleus_integral(orbital1.get_primitive(i), orbital2.get_primitive(j), nucleus_position);
        }
    }

    return integral * orbital1.normalization() * orbital2.normalization();
}