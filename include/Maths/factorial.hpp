#pragma once

long int factorial(const int n) {
    if (n == 0) { return 1; }

    long int result = 1;

    for (long int i = 1; i <= n; i++) {
        result *= i;
    }

    return result;
}