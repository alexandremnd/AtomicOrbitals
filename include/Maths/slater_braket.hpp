#pragma once

#include <boost/math/special_functions/hypergeometric_pFq.hpp>
#include <cassert>
#include <cmath>

#include "Orbitals/slater_primitive.hpp"
#include "Maths/wigner_3j.hpp"

/**
 * @brief Computes the radial integral of the product of four Slater primitives
 *
 * Integral is defined as:
    \f[
        I_1(L) = \int_\mathbb{R} \int_\mathbb{R}
            dr_1 dr_2
            R_i^* (r_1) R_j^* (r_2) R_k (r_1) R_l (r_2)
            \frac{r_<^L}{r_>^{L+1}}
            r_1^2 r_2^2
    \f]
 * where \f$ r_< = \min(r_1, r_2) \f$ and \f$ r_> = \max(r_1, r_2) \f$
 *
 * Reference: https://en.wikipedia.org/wiki/Laplace_expansion_(potential)
 *
 * @param orbital_i First Slater primitive
 * @param orbital_k Third Slater primitive sharing the coordinate with orbital_i
 * @param orbital_j Second Slater primitive
 * @param orbital_l Fourth Slater primitive sharing the coordinate with
 * orbital_j
 * @param L Legendre polynomial order of coulomb repulsion (1/r is expanded as
 Σ_L P_L(cos γ))
 * @return double Evaluation of the radial integral
 */
inline double radial_integral(const SlaterPrimitive &orbital_i,
                              const SlaterPrimitive &orbital_k,
                              const SlaterPrimitive &orbital_j,
                              const SlaterPrimitive &orbital_l, int L) {
    assert(L >= 0);
    assert(L <= std::min({orbital_i.l() + orbital_k.l(),
                          orbital_j.l() + orbital_l.l()}));

    double n_ik = orbital_i.n() + orbital_k.n();
    double n_jl = orbital_j.n() + orbital_l.n();

    double alpha_ik = orbital_i.alpha() + orbital_k.alpha();
    double alpha_jl = orbital_j.alpha() + orbital_l.alpha();

    double prefactor = boost::math::tgamma(n_ik + n_jl + 1) *
                       std::pow(alpha_ik + alpha_jl, -n_ik - n_jl - 1);

    double term = prefactor *
                  boost::math::hypergeometric_pFq<double>(
                      {1, n_ik + n_jl + 1}, {n_jl + L + 2},
                      alpha_jl / (alpha_ik + alpha_jl)) /
                  (n_jl + L + 1);
    term += prefactor *
            boost::math::hypergeometric_pFq<double>(
                {1, n_ik + n_jl + 1}, {n_ik + L + 2},
                alpha_ik / (alpha_ik + alpha_jl)) /
            (n_ik + L + 1);

    return term;
}

/**
 * @brief Computes the angular integral of the product of two Slater primitives
 *
 * Integral is defined as:
 *  \f[
 *       I_2(L, M) = \int d\Omega
 *       Y_{l_i}^{m_i} (\theta_1, \phi_1)^* Y_{l_k}^{m_k} (\theta_1, \phi_1)
 * Y_{L}^{M} (\theta_1, \phi_1)
 *  \f]
 * where \f$ Y_{l_1}^{m_1} (\theta_1, \phi_1) \f$ and \f$ Y_{l_2}^{m_2}
 * (\theta_1, \phi_1) \f$ are the spherical harmonics of the first and second
 * Slater primitives, respectively
 *
 * Reference: https://en.wikipedia.org/wiki/Laplace_expansion_(potential)
 *
 * @param orbital1 First Slater primitive
 * @param orbital2 Second Slater primitive sharing the coordinate with orbital1
 * @param L Legendre polynomial order of coulomb repulsion (1/r is expanded as
 * Σ_L P_L(cos γ)) (0 <= L)
 * @param M Expansion order of a Legendre polynomial of order L on the spherical
 * harmonics. (|M| <= L
 * @return double Evaluation of the angular integral
 */
inline double angular_integral(const SlaterPrimitive &orbital_i,
                               const SlaterPrimitive &orbital_k, int L, int M) {
    assert(L >= 0);
    assert(std::abs(M) <= L);

    double phase_factor = (orbital_i.m() % 2 == 0) ? 1. : -1.;

    return phase_factor *
           std::sqrt((2 * orbital_i.l() + 1) * (2 * orbital_k.l() + 1) *
                     (2 * L + 1) / (4 * M_PI)) *
           Math::wigner_3j(orbital_i.l(), orbital_k.l(), L, 0, 0, 0) *
           Math::wigner_3j(orbital_i.l(), orbital_k.l(), L, -orbital_i.m(),
                           orbital_k.m(), M);
}