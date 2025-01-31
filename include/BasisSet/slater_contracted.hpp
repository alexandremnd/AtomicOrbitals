#pragma once

#include <vector>

#include "BasisSet/orbital.hpp"
#include "BasisSet/slater_primitive.hpp"

class ContractedSlater final : public Orbital {
public:
    ContractedSlater() = default;
    ContractedSlater(size_t size) : m_coefficients(size), m_primitives(size) {}
    ContractedSlater(std::vector<double> &coefficients, std::vector<SlaterPrimitive> &primitives) : m_coefficients(coefficients), m_primitives(primitives) {
        update_normalization();
    }

    ContractedSlater(const ContractedSlater &other) : m_coefficients(other.m_coefficients), m_primitives(other.m_primitives) {
        m_normalization_constant = other.m_normalization_constant;
    }

    ContractedSlater(const ContractedSlater &&other) : m_coefficients(std::move(other.m_coefficients)), m_primitives(std::move(other.m_primitives)) {
        m_normalization_constant = other.m_normalization_constant;
    }


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