#pragma once

#include <vector>

#include "BasisSet/contracted_orbital.interface.hpp"
#include "BasisSet/slater_primitive.hpp"
#include "BasisSet/gaussian_primitive.hpp"
#include "concepts.hpp"


template <DerivedFromOrbital PrimitiveType>
ContractedOrbital<PrimitiveType>::ContractedOrbital(std::vector<double> &coefficients, std::vector<PrimitiveType> &primitives) : m_coefficients(coefficients), m_primitives(primitives) {
    if (coefficients.size() != primitives.size()) {
        throw std::length_error("ContractedOrbital: The number of coefficients and primitives must be equal.");
    }

    update_normalization();
}

template <DerivedFromOrbital PrimitiveType>
ContractedOrbital<PrimitiveType>::ContractedOrbital(const ContractedOrbital<PrimitiveType> &other) : m_coefficients(other.m_coefficients), m_primitives(other.m_primitives) {
    m_normalization_constant = other.m_normalization_constant;
}

template <DerivedFromOrbital PrimitiveType>
ContractedOrbital<PrimitiveType>::ContractedOrbital(const ContractedOrbital<PrimitiveType> &&other) : m_coefficients(std::move(other.m_coefficients)), m_primitives(std::move(other.m_primitives)) {
    m_normalization_constant = other.m_normalization_constant;
}

template <DerivedFromOrbital PrimitiveType>
void ContractedOrbital<PrimitiveType>::reserve(std::size_t size) {
    m_coefficients.reserve(size);
    m_primitives.reserve(size);
}

template <DerivedFromOrbital PrimitiveType>
void ContractedOrbital<PrimitiveType>::add_primitive(double coefficient, const PrimitiveType &primitive) {
    m_coefficients.push_back(coefficient);
    m_primitives.push_back(primitive);
    update_normalization();
}

template <DerivedFromOrbital PrimitiveType>
void ContractedOrbital<PrimitiveType>::update_normalization() {
    m_normalization_constant = 0.0;

    for (size_t i = 0; i < m_primitives.size(); i++) {
        for (size_t j = i; j < m_primitives.size(); j++) {
            if (i == j) {
                m_normalization_constant += m_coefficients[i] * m_coefficients[j] * overlap_integral(m_primitives[i], m_primitives[j]);
                continue;
            }

            m_normalization_constant += 2 * m_coefficients[i] * m_coefficients[j] * overlap_integral(m_primitives[i], m_primitives[j]);
        }
    }

    m_normalization_constant = 1.0 / std::sqrt(m_normalization_constant);
}

template <DerivedFromOrbital PrimitiveType>
inline size_t ContractedOrbital<PrimitiveType>::size() const {
    return m_primitives.size();
}

template <DerivedFromOrbital PrimitiveType>
inline const PrimitiveType& ContractedOrbital<PrimitiveType>::get_primitive(int i) const {
    return m_primitives[i];
}

template <DerivedFromOrbital PrimitiveType>
inline double ContractedOrbital<PrimitiveType>::get_coefficient(int i) const {
    return m_coefficients[i];
}


// ===============================================================================================
// ============================ Matrix element for Contracted Orbital ============================
// ===============================================================================================

template <DerivedFromOrbital PrimitiveType>
double overlap_integral(const ContractedOrbital<PrimitiveType>& orbital1, const ContractedOrbital<PrimitiveType>& orbital2) {
    double integral = 0.0;

    for (size_t i = 0; i < orbital1.size(); i++) {
        for (size_t j = 0; j < orbital2.size(); j++) {
            integral += orbital1.get_coefficient(i) * orbital2.get_coefficient(j) *
                        overlap_integral(orbital1.get_primitive(i), orbital2.get_primitive(j));
        }
    }

    return integral * orbital1.normalization() * orbital2.normalization();
}

template <DerivedFromOrbital PrimitiveType>
double laplacian_integral(const ContractedOrbital<PrimitiveType>& orbital1,
                        const ContractedOrbital<PrimitiveType>& orbital2) {
    double integral = 0.0;

    for (size_t i = 0; i < orbital1.size(); i++) {
        for (size_t j = 0; j < orbital2.size(); j++) {
            integral += orbital1.get_coefficient(i) * orbital2.get_coefficient(j) *
                        laplacian_integral(orbital1.get_primitive(i), orbital2.get_primitive(j));
        }
    }

    return integral * orbital1.normalization() * orbital2.normalization();
}

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

template <DerivedFromOrbital PrimitiveType>
double electron_electron_integral(const ContractedOrbital<PrimitiveType>& orbital1,
                                const ContractedOrbital<PrimitiveType>& orbital2,
                                const ContractedOrbital<PrimitiveType>& orbital3,
                                const ContractedOrbital<PrimitiveType>& orbital4) {
    double integral = 0.0;

    for (size_t i = 0; i < orbital1.size(); i++) {
        for (size_t j = 0; j < orbital2.size(); j++) {
            for (size_t k = 0; k < orbital3.size(); k++) {
                for (size_t l = 0; l < orbital4.size(); l++) {
                    integral += orbital1.get_coefficient(i) * orbital2.get_coefficient(j) * orbital3.get_coefficient(k) * orbital4.get_coefficient(l) *
                                electron_electron_integral(orbital1.get_primitive(i), orbital2.get_primitive(j), orbital3.get_primitive(k), orbital4.get_primitive(l));
                }
            }
        }
    }

    return integral * orbital1.normalization() * orbital2.normalization() * orbital3.normalization() * orbital4.normalization();
}