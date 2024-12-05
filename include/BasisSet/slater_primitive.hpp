#pragma once

#include <boost/math/special_functions/gamma.hpp>
#include <iostream>
#include <stdexcept>
#include <cassert>

class SlaterPrimitive {
public:
    SlaterPrimitive(const int n, const int l, const int m, const double alpha) : m_alpha(alpha), m_n(n), m_l(l), m_m(m) {
        assert(n > 0);
        assert(0 <= l || l < n);
        assert(abs(m) <= l);
        normalization_constant = std::pow(2 * alpha, n + 0.5) / std::sqrt(boost::math::tgamma(2 * n + 1));
    }

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