#pragma once

#include <vector>

#include "concepts.hpp"

template <DerivedFromOrbital PrimitiveType>
class ContractedOrbital final : public Orbital {
public:
    ContractedOrbital() = default;

    /**
     * @param size Number of primitives for the linear combination.
     */
    ContractedOrbital(size_t size) : m_coefficients(size), m_primitives(size) {}

    /**
     * @param coefficients Weights for each primitive in the linear combination.
     * @param primitives Primitive basis functions.
     * @throws std::length_error if the number of coefficients and primitives are not equal.
     */
    ContractedOrbital(std::vector<double> &coefficients, std::vector<ContractedOrbital<PrimitiveType>> &primitives) : m_coefficients(coefficients), m_primitives(primitives) {
        if (coefficients.size() != primitives.size()) {
            throw std::length_error("ContractedOrbital: The number of coefficients and primitives must be equal.");
        }

        update_normalization();
    }

    ContractedOrbital(const ContractedOrbital<PrimitiveType> &other) : m_coefficients(other.m_coefficients), m_primitives(other.m_primitives) {
        m_normalization_constant = other.m_normalization_constant;
    }

    ContractedOrbital(const ContractedOrbital<PrimitiveType> &&other) : m_coefficients(std::move(other.m_coefficients)), m_primitives(std::move(other.m_primitives)) {
        m_normalization_constant = other.m_normalization_constant;
    }

    /**
     * @brief Allocate memory for primitives and coefficients.
     *
     * @param size Number of primitives for the linear combination.
     */
    void reserve(std::size_t size) {
        m_coefficients.reserve(size);
        m_primitives.reserve(size);
    }

    void add_primitive(double coefficient, const PrimitiveType &primitive) {
        m_coefficients.push_back(coefficient);
        m_primitives.push_back(primitive);
        update_normalization();
    }

    template <typename... Args>
    void add_primitive(double coefficient, Args&&... args) {
        m_coefficients.push_back(coefficient);
        m_primitives.emplace_back(std::forward<Args>(args)...);
        update_normalization();
    }


    void update_normalization() {
        m_normalization_constant = 0.0;

        for (size_t i = 0; i < m_primitives.size(); i++) {
            for (size_t j = i; j < m_primitives.size(); j++) {
                m_normalization_constant += 2 * m_coefficients[i] * m_coefficients[j] * overlap_integral(m_primitives[i], m_primitives[j]);
            }
        }

        m_normalization_constant = 1.0 / std::sqrt(m_normalization_constant);
    }

    inline const PrimitiveType& get_primitive(int i) const {
        return m_primitives[i];
    }

    inline double get_coefficient(int i) const {
        return m_coefficients[i];
    }

    void set_position(Eigen::Vector3d position) override {
        this->Orbital::set_position(position);

        for (PrimitiveType &primitive : m_primitives) {
            primitive.set_position(position);
        }
    }

private:
    std::vector<double> m_coefficients;
    std::vector<PrimitiveType> m_primitives;
};

DECLARE_EXTERN_PRIMITIVE(ContractedOrbital)

typedef ContractedOrbital<GaussianPrimitive> ContractedGaussian;
typedef ContractedOrbital<SlaterPrimitive> ContractedSlater;