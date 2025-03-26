#include <vector>
#include "Orbitals/contracted_gaussian.hpp"
#include "Orbitals/gaussian_primitive.hpp"

ContractedGaussian::ContractedGaussian(
    std::vector<double> &coefficients,
    std::vector<GaussianPrimitive> &primitives)
    : m_coefficients(coefficients), m_primitives(primitives) {
    if (coefficients.size() != primitives.size()) {
        throw std::length_error(
            "ContractedGaussian: The number of coefficients "
            "and primitives must be equal.");
    }

    update_normalization();
}

ContractedGaussian::ContractedGaussian(const ContractedGaussian &other)
    : m_coefficients(other.m_coefficients), m_primitives(other.m_primitives) {
    m_constant = other.m_constant;
}

ContractedGaussian::ContractedGaussian(const ContractedGaussian &&other)
    : m_coefficients(std::move(other.m_coefficients)),
      m_primitives(std::move(other.m_primitives)) {
    m_constant = other.m_constant;
}

void ContractedGaussian::reserve(std::size_t size) {
    m_coefficients.reserve(size);
    m_primitives.reserve(size);
}

void ContractedGaussian::add_primitive(double coefficient,
                                       const GaussianPrimitive &primitive) {
    m_coefficients.push_back(coefficient);
    m_primitives.push_back(primitive);
}

void ContractedGaussian::update_normalization() {
    m_constant = 0.0;

    for (size_t i = 0; i < m_primitives.size(); i++) {
        for (size_t j = i; j < m_primitives.size(); j++) {
            if (i == j) {
                m_constant +=
                    m_coefficients[i] * m_coefficients[j] *
                    overlap_integral(m_primitives[i], m_primitives[j]);
                continue;
            }

            m_constant += 2 * m_coefficients[i] * m_coefficients[j] *
                          overlap_integral(m_primitives[i], m_primitives[j]);
        }
    }

    m_constant = 1.0 / std::sqrt(m_constant);
}

size_t ContractedGaussian::size() const { return m_primitives.size(); }

const GaussianPrimitive &ContractedGaussian::get_primitive(int i) const {
    return m_primitives[i];
}

double ContractedGaussian::get_coefficient(int i) const {
    return m_coefficients[i];
}

// ===============================================================================================
//                         Matrix element for Contracted Orbital
// ===============================================================================================

double overlap_integral(const ContractedGaussian &orbital1,
                        const ContractedGaussian &orbital2) {
    double integral = 0.0;

    for (size_t i = 0; i < orbital1.size(); i++) {
        for (size_t j = 0; j < orbital2.size(); j++) {
            double prefactor = orbital1.get_primitive(i).constant() *
                               orbital2.get_primitive(j).constant();

            integral += orbital1.get_coefficient(i) *
                        orbital2.get_coefficient(j) * prefactor *
                        overlap_integral(orbital1.get_primitive(i),
                                         orbital2.get_primitive(j));
        }
    }

    return integral;
}

double laplacian_integral(const ContractedGaussian &orbital1,
                          const ContractedGaussian &orbital2) {
    double integral = 0.0;

    for (size_t i = 0; i < orbital1.size(); i++) {
        for (size_t j = 0; j < orbital2.size(); j++) {
            double prefactor = orbital1.get_primitive(i).constant() *
                               orbital2.get_primitive(j).constant();

            integral += orbital1.get_coefficient(i) *
                        orbital2.get_coefficient(j) * prefactor *
                        laplacian_integral(orbital1.get_primitive(i),
                                           orbital2.get_primitive(j));
        }
    }

    return integral;
}

double electron_nucleus_integral(const ContractedGaussian &orbital1,
                                 const ContractedGaussian &orbital2,
                                 const Eigen::Vector3d &nucleus_position) {
    double integral = 0.0;
    for (size_t i = 0; i < orbital1.size(); i++) {
        for (size_t j = 0; j < orbital2.size(); j++) {
            double prefactor = orbital1.get_primitive(i).constant() *
                               orbital2.get_primitive(j).constant();

            integral += orbital1.get_coefficient(i) *
                        orbital2.get_coefficient(j) * prefactor *
                        electron_nucleus_integral(orbital1.get_primitive(i),
                                                  orbital2.get_primitive(j),
                                                  nucleus_position);
        }
    }

    return integral;
}

double electron_electron_integral(const ContractedGaussian &orbital1,
                                  const ContractedGaussian &orbital2,
                                  const ContractedGaussian &orbital3,
                                  const ContractedGaussian &orbital4) {
    double integral = 0.0;

    for (size_t i = 0; i < orbital1.size(); i++) {
        for (size_t j = 0; j < orbital2.size(); j++) {
            for (size_t k = 0; k < orbital3.size(); k++) {
                for (size_t l = 0; l < orbital4.size(); l++) {
                    double prefactor = orbital1.get_primitive(i).constant() *
                                       orbital2.get_primitive(j).constant() *
                                       orbital3.get_primitive(k).constant() *
                                       orbital4.get_primitive(l).constant();

                    integral +=
                        orbital1.get_coefficient(i) *
                        orbital2.get_coefficient(j) *
                        orbital3.get_coefficient(k) *
                        orbital4.get_coefficient(l) * prefactor *
                        electron_electron_integral(orbital1.get_primitive(i),
                                                   orbital2.get_primitive(j),
                                                   orbital3.get_primitive(k),
                                                   orbital4.get_primitive(l));
                }
            }
        }
    }

    return integral;
}