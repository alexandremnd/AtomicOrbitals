#pragma once

#include <vector>
#include "BasisSet/gaussian_primitive.hpp"

class ContractedGaussian {
public:
    ContractedGaussian(const std::vector<double>& coefficients,
                       const std::vector<GaussianPrimitive>& primitives)
        : m_coefficients(coefficients), m_primitives(primitives) {}

    double evaluate(double x, double y, double z) const;

    int get_gaussian_count() const {
        return m_primitives.size();
    }

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