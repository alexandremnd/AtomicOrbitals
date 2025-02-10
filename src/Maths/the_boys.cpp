#include <boost/math/special_functions/gamma.hpp>
#include <cmath>
#include <vector>

#include "Maths/the_boys.hpp"

double TheBoys_Function(int &n, double &x) {
    if (x <= 1e-3) {
        return 1. / (2 * n + 1);
    }
    return boost::math::tgamma_lower(n + (1. / 2), x) /
           (2. *
            std::pow(x, n + (1. / 2))); // Return the values of The Boys
                                        // function of order n evaluated at x
}

std::vector<double> TheBoys_Recurrence(int &n, double &x) {
    std::vector<double> Fn_values(n + 1);
    for (int m = n; m >= 0; m--) {
        if (m == n) {
            Fn_values[m] = TheBoys_Function(m, x);
        }

        else {
            Fn_values[m] =
                (2 * x * Fn_values[m + 1] + std::exp(-x)) / (2. * m + 1.);
        }
    }
    return Fn_values; // Return the values of Fm(x) for 0<=m<=n like this
                      // [F0,F1,...,Fn]
}