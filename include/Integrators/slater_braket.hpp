#pragma once

#include <boost/math/special_functions/hypergeometric_pFq.hpp>
#include <cmath>

#include "BasisSet/slater_primitive.hpp"
#include "Maths/wigner_3j.hpp"

/**
 * @brief Compute the radial integral of the product of four Slater primitives
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
 * @param orbital1 First Slater primitive
 * @param orbital2 Second Slater primitive
 * @param orbital3 Third Slater primitive
 * @param orbital4 Fourth Slater primitive
 * @param L Angular momentum quantum number
 * @return double Evaluation of the radial integral
 */
inline double radial_integral(const SlaterPrimitive &orbital1, const SlaterPrimitive &orbital2, const SlaterPrimitive &orbital3, const SlaterPrimitive &orbital4, const int L) {
    double n_ik = orbital1.n() + orbital3.n();
    double n_jl = orbital2.n() + orbital4.n();

    double alpha_ik = orbital1.alpha() + orbital3.alpha();
    double alpha_jl = orbital2.alpha() + orbital4.alpha();

    double prefactor = boost::math::tgamma(n_ik + n_jl + 1) * std::pow(alpha_ik + alpha_jl, -n_ik - n_jl - 1);

    double term = prefactor * boost::math::hypergeometric_pFq<double>({1, n_ik + n_jl + 1}, {n_jl + L + 2}, alpha_jl / (alpha_ik + alpha_jl)) / (n_jl + L + 1);
    term += prefactor * boost::math::hypergeometric_pFq<double>({1, n_ik + n_jl + 1}, {n_ik + L + 2}, alpha_ik / (alpha_ik + alpha_jl)) / (n_ik + L + 1);

    return term;
}

/**
 * @brief Compute the angular integral of the product of two Slater primitives
 *
 * Integral is defined as:
 *  \f[
 *       I_2(L, M) = \int d\Omega
 *       Y_{l_1}^{m_1} (\theta_1, \phi_1)^* Y_{l_2}^{m_2} (\theta_1, \phi_1) Y_{L}^{M} (\theta_1, \phi_1)
 *  \f]
 * where \f$ Y_{l_1}^{m_1} (\theta_1, \phi_1) \f$ and \f$ Y_{l_2}^{m_2} (\theta_1, \phi_1) \f$ are the spherical harmonics of the first and second Slater primitives, respectively
 *
 * @param orbital1 First Slater primitive
 * @param orbital2 Second Slater primitive
 * @param L Angular momentum quantum number
 * @param M Magnetic quantum number
 * @return double Evaluation of the angular integral
 */
inline double angular_integral(const SlaterPrimitive &orbital1, const SlaterPrimitive &orbital2, const int L, const int M) {
    double phase_factor = (orbital1.m() % 2 == 0) ? 1. : -1.;

    return phase_factor * std::sqrt( (2 * orbital1.l() + 1) * (2 * orbital2.l() + 1) * (2 * L + 1) / (4 * M_PI) )
    * Math::wigner_3j(orbital1.l(), orbital2.l(), L, 0, 0, 0)
    * Math::wigner_3j(orbital1.l(), orbital2.l(), L, -orbital1.m(), orbital2.m(), M);
}