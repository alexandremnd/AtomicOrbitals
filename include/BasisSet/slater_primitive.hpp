#pragma once

#include <boost/math/special_functions/gamma.hpp>
#include <stdexcept>
#include <cassert>

#include "BasisSet/orbital.hpp"


/**
 * @brief Slater type orbital (STO) representation

    It is the braket representation \f$ | R_n^\alpha l m \rangle \f$ where :
    \f[
        R_n^\alpha (r) = r^{n-1}  e^{-\alpha r}
    \f]
    and \f$ | l m \rangle \f$ are angular part of the wavefunction which are the spherical harmonics.
 **/
class SlaterPrimitive : public Orbital {
public:
    SlaterPrimitive() = default;
    ~SlaterPrimitive() = default;

    /**
     * Throws an exception (std::invalid_argument) if one or more parameters are invalid.
     *
     * Note: The normalization constant only normalize the radial part of the wavefunction as
     * we expressed braket STO without normalization constant.
     * The angular part is normalized in the spherical harmonics.
     *
     * @param n Principal quantum number (n > 0)
     * @param l Secondary quantum number (0 < l < n)
     * @param m Magnetic quantum number (-l <= m <= l)
     * @param alpha Exponential decay constant (alpha > 0)
     */
    SlaterPrimitive(const int n, const int l, const int m, const double alpha) : m_n(n), m_l(l), m_m(m), m_alpha(alpha) {
        check_parameters(n, l, m, alpha);
        update_normalization_constant();
    }

    /**
     * Throws an exception (std::invalid_argument) if the parameter is invalid.
     *
     * Updates the normalization constant if parameter is accepted.
     *
     * @param alpha Exponential decay constant (alpha > 0)
     */
    void set_alpha(const double alpha) {
        check_parameters(m_n, m_l, m_m, alpha);

        m_alpha = alpha;
        update_normalization_constant();
    }

    /**
     * Throws an exception (std::invalid_argument) if the parameter is invalid.
     *
     * Updates the normalization constant if parameter is accepted.
     *
     * @param n Principal quantum number (n > 0)
     */
    void set_n(const int n) {
        check_parameters(n, m_l, m_m, m_alpha);

        m_n = n;
        update_normalization_constant();
    }

    /**
     * Throws an exception (std::invalid_argument) if the parameter is invalid.
     *
     * @param l Azimuthal quantum number (0 < l < n)
     */
    void set_l(const int l) {
        check_parameters(m_n, l, m_m, m_alpha);
        m_l = l;
    }

    /**
     * Throws an exception (std::invalid_argument) if the parameter is invalid.
     *
     * @param m Magnetic quantum number (-l <= m <= l)
     */
    void set_m(const int m) {
        check_parameters(m_n, m_l, m, m_alpha);
        m_m = m;
    }

    inline int n() const { return m_n; }
    inline int l() const { return m_l; }
    inline int m() const { return m_m; }
    inline double alpha() const { return m_alpha; }

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

    void update_normalization_constant() {
        m_normalization_constant = std::pow(2 * m_alpha, m_n + 0.5) / std::sqrt(boost::math::tgamma(2 * m_n + 1));
    }

    int m_n, m_l, m_m;
    double m_alpha;
};