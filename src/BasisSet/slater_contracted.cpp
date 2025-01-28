#include "BasisSet/slater_contracted.hpp"
#include "Integrators/overlap_integral.hpp"

void ContractedSlater::reserve(std::size_t size) {
    m_coefficients.reserve(size);
    m_primitives.reserve(size);
}

void ContractedSlater::add_primitive(const double coefficient, const SlaterPrimitive &primitive) {
    m_coefficients.push_back(coefficient);
    m_primitives.push_back(primitive);
}

void ContractedSlater::add_primitive(double coefficient, const int n, const int l, const int m, const double alpha) {
    m_coefficients.push_back(coefficient);
    m_primitives.emplace_back(n, l, m, alpha);
}

void ContractedSlater::update_normalization() {
    m_normalization_constant = 0.0;

    for (size_t i = 0; i < m_primitives.size(); i++) {
        for (size_t j = i; j < m_primitives.size(); j++) {
            m_normalization_constant += 2 * m_coefficients[i] * m_coefficients[j] * overlap_integral(m_primitives[i], m_primitives[j]);
        }
    }

    m_normalization_constant = 1.0 / std::sqrt(this->m_normalization_constant);
}