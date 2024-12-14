#pragma once

#include <boost/math/special_functions/gamma.hpp>
#include <cmath>
#include <algorithm>

namespace Math
{
    /**
     * @brief Computes the Wigner's 3-j symbol
     *
     * Definition given in: https://en.wikipedia.org/wiki/3-j_symbol
     *
     * @param j1 First angular momentum quantum number
     * @param j2 Second angular momentum quantum number
     * @param j3 Third angular momentum quantum number
     * @param m1 Projection of the first angular momentum quantum number
     * @param m2 Projection of the second angular momentum quantum number
     * @param m3 Projection of the third angular momentum quantum number
     * @return double The Wigner's 3-j symbol value
     */
    inline double wigner_3j(const int j1, const int j2, const int j3, const int m1, const int m2, const int m3) {
        if (m1 + m2 + m3 != 0) {
            return 0.0;
        }

        if (j1 + j2 - j3 < 0 || j1 - j2 + j3 < 0 || -j1 + j2 + j3 < 0) {
            return 0.0;
        }

        if (j1 - m1 < 0 || j2 - m2 < 0 || j3 - m3 < 0 ) {
            return 0.0;
        }

        double phase_factor = ((j1 - j2 - m3) % 2) == 0 ? 1. : -1.;

        // Summation bounds for sum from k=K to N
        long int N = std::min({j1 + j2 - j3, j1 - m1, j2 + m2});
        long int K = std::max({0, j2 - j3 - m1, j1 - j3 + m2});

        double t1 = j1 + j2 - j3;
        double t2 = j1 - m1;
        double t3 = j2 + m2;
        double t4 = j3 - j2 + m1;
        double t5 = j3 - j1 - m2;

        double sum = ((K % 2) == 0) ? 1.0 : -1.0;
        double term = 1.
                    / boost::math::tgamma(K + 1)
                    / boost::math::tgamma(t1 - K + 1)
                    / boost::math::tgamma(t2 - K + 1)
                    / boost::math::tgamma(t3 - K + 1)
                    / boost::math::tgamma(t4 + K + 1)
                    / boost::math::tgamma(t5 + K + 1);
        sum *= term;

        // Recurrence relation used on the sum over L in 3j symbol
        for (int i = K; i < N; i++) {
            term *= -1.0/(i + 1);
            term *= (t1 - i) * (t2 - i) * (t3 - i);
            term /= (t4 + i + 1) * (t5 + i + 1);
            sum += term;
        }

        double prefactor = phase_factor * std::sqrt(
            boost::math::tgamma(j1 + j2 - j3 + 1)
            * boost::math::tgamma(j1 - j2 + j3 + 1)
            * boost::math::tgamma(-j1 + j2 + j3 + 1)
            / boost::math::tgamma(j1 + j2 + j3 + 2)
            * boost::math::tgamma(j1 - m1 + 1)
            * boost::math::tgamma(j1 + m1 + 1)
            * boost::math::tgamma(j2 - m2 + 1)
            * boost::math::tgamma(j2 + m2 + 1)
            * boost::math::tgamma(j3 - m3 + 1)
            * boost::math::tgamma(j3 + m3 + 1)
        );

        return prefactor * sum;
    }
}