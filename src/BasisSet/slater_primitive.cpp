#pragma once

#include <boost/math/special_functions/gamma.hpp>
#include <stdexcept>

#include "BasisSet/slater_primitive.hpp"

SlaterPrimitive::SlaterPrimitive(const int n, const int l, const int m, const double alpha) : m_n(n), m_l(l), m_m(m), m_alpha(alpha) {
    check_parameters(n, l, m, alpha);
    update_normalization_constant();
}

void SlaterPrimitive::set_alpha(const double alpha) {
    check_parameters(m_n, m_l, m_m, alpha);

    m_alpha = alpha;
    update_normalization_constant();
}


void SlaterPrimitive::set_n(const int n) {
    check_parameters(n, m_l, m_m, m_alpha);

    m_n = n;
    update_normalization_constant();
}

/**
    * Throws an exception (std::invalid_argument) if the parameter is invalid.
    *
    * @param l Azimuthal quantum number (0 < l < n)
    */
void SlaterPrimitive::set_l(const int l) {
    check_parameters(m_n, l, m_m, m_alpha);
    m_l = l;
}

/**
    * Throws an exception (std::invalid_argument) if the parameter is invalid.
    *
    * @param m Magnetic quantum number (-l <= m <= l)
    */
void SlaterPrimitive::set_m(const int m) {
    check_parameters(m_n, m_l, m, m_alpha);
    m_m = m;
}

void SlaterPrimitive::check_parameters(const int n, const int l, const int m, const double alpha) {
    if (n < 1) {
        throw std::invalid_argument("SlaterPrimitive: n must be a positive integer greater than 0");
    }
    if (l < 0 || l >= n) {
        throw std::invalid_argument("SlaterPrimitive: l must be between 0 and n - 1");
    }
    if (abs(m) > l) {
        throw std::invalid_argument("SlaterPrimitive: m must be between -l and l");
    }
    if (alpha <= 0) {
        throw std::invalid_argument("SlaterPrimitive: alpha must be a positive real number");
    }
}

void SlaterPrimitive::update_normalization_constant() {
    m_normalization_constant = std::pow(2 * m_alpha, m_n + 0.5) / std::sqrt(boost::math::tgamma(2 * m_n + 1));
}