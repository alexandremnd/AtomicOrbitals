#pragma once

#include <vector>
#include "include/BasisSet/gaussian_primitive.hpp"

class ContractedGaussian {
public:
    ContractedGaussian(const std::vector<double>& coefficients,
                       const std::vector<GaussianPrimitive>& primitives)
        : m_coefficients(coefficients), m_primitives(primitives) {}

private:
    std::vector<double> m_coefficients;
    std::vector<GaussianPrimitive> m_primitives;

}