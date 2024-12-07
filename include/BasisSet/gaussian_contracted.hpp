#pragma once

#include <vector>
#include "include/BasisSet/gaussian_primitive.hpp"

class ContractedGaussian {
public:
    ContractedGaussian(const std::vector<double>& coefficients,
                       const std::vector<GaussianPrimitive>& primitives)
        : m_coefficients(coefficients), m_primitives(primitives) {}

    double evaluate(double x, double y, double z) const {}

    double get_coefficient(int i) const {
        return m_coefficients[i];
    }

    const GaussianPrimitive& get_primitive(int i) const {
        return m_primitives[i];
    }

private:
    std::vector<double> m_coefficients;
    std::vector<GaussianPrimitive> m_primitives;
};