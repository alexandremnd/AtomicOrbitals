#pragma once

#include <boost/math/special_functions/hypergeometric_pFq.hpp>
#include <cmath>

#include "include/BasisSet/slater_primitive.hpp"
#include "include/Maths/wigner_3j.hpp"

inline double radial_integral(const SlaterPrimitive &orbital1, const SlaterPrimitive &orbital2, const SlaterPrimitive &orbital3, const SlaterPrimitive &orbital4, const int L) {
    double n_ik = orbital1.n() + orbital3.n();
    double n_jl = orbital2.n() + orbital4.n();

    double alpha_ik = orbital1.alpha() + orbital3.alpha();
    double alpha_jl = orbital2.alpha() + orbital4.alpha();

    double prefactor = orbital1.normalization() * orbital2.normalization() * orbital3.normalization() * orbital4.normalization() * boost::math::tgamma(n_ik + n_jl + 1) * std::pow(alpha_ik + alpha_jl, -n_ik - n_jl - 1);

    double term = prefactor * boost::math::hypergeometric_pFq<double>({1, n_ik + n_jl + 1}, {n_jl + L + 2}, alpha_jl / (alpha_ik + alpha_jl)) / (n_jl + L + 1);
    term += prefactor * boost::math::hypergeometric_pFq<double>({1, n_ik + n_jl + 1}, {n_ik + L + 2}, alpha_ik / (alpha_ik + alpha_jl)) / (n_ik + L + 1);

    return term;
}

inline double angular_integral(const SlaterPrimitive &orbital1, const SlaterPrimitive &orbital2, const int L, const int M) {
    return std::sqrt( (2 * orbital1.l() + 1) * (2 * orbital2.l()) * (2 * L + 1) / (4 * M_PI) ) * Math::wigner_3j(orbital1.l(), orbital2.l(), L, 0, 0, M) * Math::wigner_3j(orbital1.l(), orbital2.l(), L, orbital1.m(), orbital2.m(), M);
}