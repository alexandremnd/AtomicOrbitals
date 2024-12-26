#pragma once

#include <boost/math/special_functions/gamma.hpp>
#include <stdexcept>
#include <cassert>


class SlaterPrimitive {
public:
    /**
     * @brief Construct a new Slater Primitive object
     *
     * Note: The normalization constant only normalize the radial part of the wavefunction as
     * we expressed braket STO without normalization constant.
     * The angular part is normalized in the spherical harmonics.
     *
     * @param n Principal quantum number
     * @param l Secondary quantum number
     * @param m Magnetic quantum number
     * @param alpha Exponential decay constant
     */
    SlaterPrimitive(const int n, const int l, const int m, const double alpha) : m_n(n), m_l(l), m_m(m), m_alpha(alpha) {
        check_parameters(n, l, m, alpha);
        compute_normalization_constant();
    }

    void set_alpha(const double alpha) {
        check_parameters(m_n, m_l, m_m, alpha);

        m_alpha = alpha;
        compute_normalization_constant();
    }

    void set_n(const int n) {
        check_parameters(n, m_l, m_m, m_alpha);

        m_n = n;
        compute_normalization_constant();
    }

    void set_l(const int l) {
        check_parameters(m_n, l, m_m, m_alpha);
        m_l = l;
    }

    void set_m(const int m) {
        check_parameters(m_n, m_l, m, m_alpha);
        m_m = m;
    }

    void increment_n() { set_n(m_n + 1); }
    void decrement_n() { set_n(m_n - 1); }

    inline int n() const { return m_n; }
    inline int l() const { return m_l; }
    inline int m() const { return m_m; }
    inline double alpha() const { return m_alpha; }
    inline double normalization() const { return m_normalization_constant; }

private:
    void check_parameters(const int n, const int l, const int m, const double alpha) {
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

    void compute_normalization_constant() {
        m_normalization_constant = std::pow(2 * m_alpha, m_n + 0.5) / std::sqrt(boost::math::tgamma(2 * m_n + 1));
    }

    int m_n, m_l, m_m;
    double m_alpha;
    double m_normalization_constant;
};