#pragma once

#include <boost/math/special_functions/gamma.hpp>
#include <iostream>
#include <stdexcept>
#include <cassert>

class SlaterPrimitive {
public:
    SlaterPrimitive(const int n, const int l, const int m, const double alpha) : m_n(n), m_l(l), m_m(m), m_alpha(alpha) {
        assert(1 <= n);
        assert(0 <= l); // We raise the constraint of l <= n - 1 to allow possible polarization
        assert(abs(m) <= l);
        normalization_constant = std::pow(2 * alpha, n + 0.5) / std::sqrt(boost::math::tgamma(2 * n + 1));
    }

    void set_n(const int n) { m_n = n; }

    void set_l(const int l) {
        if (l < 0) {
            throw std::invalid_argument("l must be a non-negative integer");
        }
        m_l = l;
    }

    void set_m(const int m) {
        if (abs(m) > m_l) {
            throw std::invalid_argument("m must be between -l and l");
        }
        m_m = m;
    }

    void increment_n() { set_n(m_n + 1); }
    void decrement_n() { set_n(m_n - 1); }

    int n() const { return m_n; }
    int l() const { return m_l; }
    int m() const { return m_m; }
    double alpha() const { return m_alpha; }
    double normalization() const { return normalization_constant; }

private:
    int m_n, m_l, m_m;
    double m_alpha;
    double normalization_constant;

};