#include <stdexcept>

#include "Integrators/electron_electron.hpp"
#include "Integrators/slater_braket.hpp"
#include "BasisSet/gaussian_primitive.hpp"
#include "BasisSet/slater_primitive.hpp"

double electron_electron_integral(const GaussianPrimitive &orbital1, const GaussianPrimitive &orbital2, const GaussianPrimitive &orbital3, const GaussianPrimitive &orbital4)
{
    throw std::logic_error("Not implemented yet");
}

double electron_electron_integral(const SlaterPrimitive &orbital1, const SlaterPrimitive &orbital2, const SlaterPrimitive &orbital3, const SlaterPrimitive &orbital4)
{
    int max_L = std::min({orbital1.l() + orbital3.l(), orbital2.l() + orbital4.l()});
    double matrix_element = 0.0;

    for (int L = 0; L <= max_L; L++) {
        double value = 0.0;

        for (int M = -L; M <= L; M++) {
            // Wigner 3-j symbol are non-zero only if the sum of the m's is zero
            if (-orbital1.m() + orbital3.m() + M != 0 || -orbital2.m() + orbital4.m() - M != 0) {
                continue;
            }

            double m_value = (M % 2 == 0) ? 1. : -1.;
            m_value *= angular_integral(orbital1, orbital3, L, M);
            m_value *= angular_integral(orbital2, orbital4, L, -M) ;

            value += m_value;
        }

        value *= 4. * M_PI / (2. * L + 1.) * radial_integral(orbital1, orbital2, orbital3, orbital4, L);
        matrix_element += value;
    }

    return matrix_element * orbital1.normalization() * orbital2.normalization() * orbital3.normalization() * orbital4.normalization();
}
