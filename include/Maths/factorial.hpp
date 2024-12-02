#pragma once

namespace Math
{
    /**
     * @brief Computes the factorial of a number
     *
     * @param n The number to compute the factorial of
     * @return long int The factorial of the number
     */
    long int factorial(const int n) {
        if (n == 0) { return 1; }

        long int result = 1;

        for (long int i = 1; i <= n; i++) {
            result *= i;
        }

        return result;
    }
}