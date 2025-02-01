#pragma once

#include <vector>
#include <stdexcept>

#include "BasisSet/orbital.hpp"
#include "BasisSet/slater_primitive.hpp"

/**
 * @brief Contracted slater type orbital representation.
    It is the linear combination of Slater primitives \f$ | R_n^\alpha l m \rangle \f$ such that :
    \f[
        \Phi_i = \sum_i c_i R_{n_i}^{\alpha_i} (r)
    \f]
    where \f$ c_i \f$ are the coefficients and \f$ R_{n_i}^{\alpha_i} (r) \f$ are the Slater primitives.
 */
class ContractedSlater final : public Orbital {
public:
    ContractedSlater() = default;

    /**
     * @param size Number of primitives for the linear combination.
     */
    ContractedSlater(size_t size) : m_coefficients(size), m_primitives(size) {}

    /**
     * @param coefficients Weights for each primitive in the linear combination.
     * @param primitives Primitive basis functions.
     * @throws std::length_error if the number of coefficients and primitives are not equal.
     */
    ContractedSlater(std::vector<double> &coefficients, std::vector<SlaterPrimitive> &primitives) : m_coefficients(coefficients), m_primitives(primitives) {
        if (coefficients.size() != primitives.size()) {
            throw std::length_error("ContractedSlater: The number of coefficients and primitives must be equal.");
        }

        update_normalization();
    }

    ContractedSlater(const ContractedSlater &other) : m_coefficients(other.m_coefficients), m_primitives(other.m_primitives) {
        m_normalization_constant = other.m_normalization_constant;
    }

    ContractedSlater(const ContractedSlater &&other) : m_coefficients(std::move(other.m_coefficients)), m_primitives(std::move(other.m_primitives)) {
        m_normalization_constant = other.m_normalization_constant;
    }

    /**
     * @brief Allocate memory for primitives and coefficients.
     *
     * @param size Number of primitives for the linear combination.
     */
    void reserve(std::size_t size);

    void add_primitive(double coefficient, const SlaterPrimitive &primitive);
    void add_primitive(double coefficient, int n, int l, int m, double alpha);
    void update_normalization();

    inline const SlaterPrimitive& get_primitive(int i) const {
        return m_primitives[i];
    }

    inline double get_coefficient(int i) const {
        return m_coefficients[i];
    }

private:
    std::vector<double> m_coefficients;
    std::vector<SlaterPrimitive> m_primitives;
};