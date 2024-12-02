#pragma once

#include "include/Maths/factorial.hpp"

#include <cmath>
#include <algorithm>

namespace Math
{
    /**
     * @brief Computes the Wigner's 3-j symbol
     *
     * Definition given in: https://en.wikipedia.org/wiki/3-j_symbol
     * Logarithmic gamma function is used to avoid (long) integer overflow, at the cost of precision.
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

        int phase_factor = ((j1 + j2 - m3 % 2) == 0) ? 1 : -1;

        // Summation bounds for sum from k=K to N
        long int N = std::min({j1 + j2 - j3, j1 - m1, j2 + m2});
        long int K = std::max({0, j2 - j3 - m1, j1 - j3 + m2});

        long int t1 = j1 + j2 - j3;
        long int t2 = j1 - m1;
        long int t3 = j2 + m2;
        long int t4 = j3 - j2 + m1;
        long int t5 = j3 - j1 - m2;

        double sum = ((K % 2) == 0) ? 1.0 : -1.0;
        double term = std::exp(
            - std::lgamma(K + 1)
            - std::lgamma(t1 - K + 1)
            - std::lgamma(t2 - K + 1)
            - std::lgamma(t3 - K + 1)
            - std::lgamma(t4 + K + 1)
            - std::lgamma(t5 + K + 1)
        );
        sum *= term;

        for (int k = K; k < N; k++) {
            term *= -1.0/(k + 1);
            term *= (t1 - k) * (t2 - k) * (t3 - k);
            term /= (t4 + k + 1) * (t5 + k + 1);
            sum += term;
        }

        double prefactor = phase_factor * std::sqrt(
            std::exp(
                std::lgamma(j1 + j2 - j3 + 1) +
                std::lgamma(j1 - j2 + j3 + 1) +
                std::lgamma(-j1 + j2 + j3 + 1) -
                std::lgamma(j1 + j2 + j3 + 2) +
                std::lgamma(j1 - m1 + 1) +
                std::lgamma(j1 + m1 + 1) +
                std::lgamma(j2 - m2 + 1) +
                std::lgamma(j2 + m2 + 1) +
                std::lgamma(j3 - m3 + 1) +
                std::lgamma(j3 + m3 + 1)
            )
        );

        return prefactor * sum;
    }
}