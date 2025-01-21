#pragma once

#include <vector>

#include "BasisSet/slater_primitive.hpp"
#include "Integrators/overlap_integral.hpp"

class ContractedSlater : public Orbital {
public:
    ContractedSlater() = default;
    ContractedSlater(std::vector<double>& coefficients, std::vector<SlaterPrimitive>& primitives) : m_coefficients(coefficients), m_primitives(primitives) {
        this->m_normalization_constant = 0.0;

        for (size_t i = 0; i < m_primitives.size(); i++) {
            for (size_t j = i; j < m_primitives.size(); j++) {
                this->m_normalization_constant += 2 * m_coefficients[i] * m_coefficients[j] * overlap_integral(m_primitives[i], m_primitives[j]);
            }
        }

        this->m_normalization_constant = 1.0 / std::sqrt(this->m_normalization_constant);
    }

private:
    std::vector<double> m_coefficients;
    std::vector<SlaterPrimitive> m_primitives;
};