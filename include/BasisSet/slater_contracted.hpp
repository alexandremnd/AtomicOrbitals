#pragma once

#include <cstddef>
#include <vector>

#include "BasisSet/orbital.hpp"
#include "BasisSet/slater_primitive.hpp"

class ContractedSlater : public Orbital {
public:
    ContractedSlater() = default;
    ContractedSlater(std::size_t size) : m_coefficients(size), m_primitives(size) {}
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
    void add_primitive(const double coefficient, const SlaterPrimitive &primitive);
    void add_primitive(double coefficient, const int n, const int l, const int m, const double alpha);
    void update_normalization();

private:
    std::vector<double> m_coefficients;
    std::vector<SlaterPrimitive> m_primitives;
};