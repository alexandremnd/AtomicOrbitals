#include <boost/math/special_functions/gamma.hpp>
#include <stdexcept>

#include "Orbitals/slater_primitive.hpp"
#include "Integrators/slater_braket.hpp"

SlaterPrimitive::SlaterPrimitive(int n, int l, int m, double alpha)
    : m_n(n), m_l(l), m_m(m), m_alpha(alpha) {
    check_parameters(n, l, m, alpha);
    update_normalization_constant();
}

void SlaterPrimitive::set_alpha(double alpha) {
    check_parameters(m_n, m_l, m_m, alpha);

    m_alpha = alpha;
    update_normalization_constant();
}

void SlaterPrimitive::set_n(int n) {
    check_parameters(n, m_l, m_m, m_alpha);

    m_n = n;
    update_normalization_constant();
}

void SlaterPrimitive::set_l(int l) {
    check_parameters(m_n, l, m_m, m_alpha);
    m_l = l;
}

void SlaterPrimitive::set_m(int m) {
    check_parameters(m_n, m_l, m, m_alpha);
    m_m = m;
}

void SlaterPrimitive::check_parameters(int n, int l, int m, double alpha) {
    if (n < 1) {
        throw std::invalid_argument(
            "SlaterPrimitive: n must be a positive integer greater than 0");
    }
    if (l < 0 || l >= n) {
        throw std::invalid_argument(
            "SlaterPrimitive: l must be between 0 and n - 1");
    }
    if (abs(m) > l) {
        throw std::invalid_argument(
            "SlaterPrimitive: m must be between -l and l");
    }
    if (alpha <= 0) {
        throw std::invalid_argument(
            "SlaterPrimitive: alpha must be a positive real number");
    }
}

void SlaterPrimitive::update_normalization_constant() {
    m_normalization_constant = std::pow(2 * m_alpha, m_n + 0.5) /
                               std::sqrt(boost::math::tgamma(2 * m_n + 1));
}

// ===============================================================================================
// ============================ Matrix element for Slater Primitive
// ==============================
// ===============================================================================================

double overlap_integral(const SlaterPrimitive &orbital1,
                        const SlaterPrimitive &orbital2, const int n_offset) {
    assert(orbital1.n() + orbital2.n() + n_offset >= 0);

    if (orbital1.l() != orbital2.l() || orbital1.m() != orbital2.m()) {
        return 0.;
    }

    return boost::math::tgamma(orbital1.n() + orbital2.n() + n_offset + 1) *
           std::pow(orbital1.alpha() + orbital2.alpha(),
                    -orbital1.n() - orbital2.n() - n_offset - 1);
}

double laplacian_integral(const SlaterPrimitive &orbital1,
                          const SlaterPrimitive &orbital2) {
    double matrix_element = orbital2.alpha() * orbital2.alpha() *
                            overlap_integral(orbital1, orbital2);
    matrix_element -= 2 * orbital2.alpha() * orbital2.n() *
                      overlap_integral(orbital1, orbital2, -1);
    matrix_element += (orbital2.n() * (orbital2.n() - 1) -
                       orbital2.l() * (orbital2.l() + 1)) *
                      overlap_integral(orbital1, orbital2, -2);

    return matrix_element;
}

double electron_nucleus_integral(const SlaterPrimitive &orbital1,
                                 const SlaterPrimitive &orbital2,
                                 const Eigen::Vector3d &nucleus_position) {
    assert(nucleus_position.isZero()); // Slater HF computation are not
                                       // implemented for multi-atoms systems

    return overlap_integral(orbital1, orbital2, -1);
}

double electron_electron_integral(const SlaterPrimitive &orbital1,
                                  const SlaterPrimitive &orbital2,
                                  const SlaterPrimitive &orbital3,
                                  const SlaterPrimitive &orbital4) {
    int max_L =
        std::min({orbital1.l() + orbital3.l(), orbital2.l() + orbital4.l()});
    double matrix_element = 0.0;

    for (int L = 0; L <= max_L; L++) {
        double value = 0.0;

        for (int M = -L; M <= L; M++) {
            // Wigner 3-j symbol are non-zero only if the sum of the m's is zero
            if (-orbital1.m() + orbital3.m() + M != 0 ||
                -orbital2.m() + orbital4.m() - M != 0) {
                continue;
            }

            double m_value = (M % 2 == 0) ? 1. : -1.;
            m_value *= angular_integral(orbital1, orbital3, L, M);
            m_value *= angular_integral(orbital2, orbital4, L, -M);

            value += m_value;
        }

        value *= 4. * M_PI / (2. * L + 1.) *
                 radial_integral(orbital1, orbital2, orbital3, orbital4, L);
        matrix_element += value;
    }

    return matrix_element;
}