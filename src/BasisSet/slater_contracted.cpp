#include "BasisSet/slater_contracted.hpp"

#include "BasisSet/slater_primitive.hpp"
#include "Integrators/overlap_integral.hpp"

void ContractedSlater::reserve(size_t size) {
    m_coefficients.reserve(size);
    m_primitives.reserve(size);
}

void ContractedSlater::add_primitive(double coefficient, const SlaterPrimitive &primitive) {
    m_coefficients.push_back(coefficient);
    m_primitives.push_back(primitive);
    update_normalization();
}

void ContractedSlater::add_primitive(double coefficient, int n, int l, int m, double alpha) {
    m_coefficients.push_back(coefficient);
    m_primitives.emplace_back(n, l, m, alpha);
    update_normalization();
}

void ContractedSlater::update_normalization() {
    m_normalization_constant = 0.0;

    for (size_t i = 0; i < m_primitives.size(); i++) {
        for (size_t j = i; j < m_primitives.size(); j++) {
            m_normalization_constant += 2 * m_coefficients[i] * m_coefficients[j] * overlap_integral(m_primitives[i], m_primitives[j]);
        }
    }

    m_normalization_constant = 1.0 / std::sqrt(m_normalization_constant);
}