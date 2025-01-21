#include "Integrators/slater_integrator.hpp"
#include "BasisSet/slater_primitive.hpp"
#include "Integrators/slater_braket.hpp"

double SlaterIntegrator::electron_electron_integral(const Orbital& o1, const Orbital& o2, const Orbital& o3, const Orbital& o4) const {
    const SlaterPrimitive& orbital1 = dynamic_cast<const SlaterPrimitive&>(o1);
    const SlaterPrimitive& orbital2 = dynamic_cast<const SlaterPrimitive&>(o2);
    const SlaterPrimitive& orbital3 = dynamic_cast<const SlaterPrimitive&>(o3);
    const SlaterPrimitive& orbital4 = dynamic_cast<const SlaterPrimitive&>(o4);

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

double SlaterIntegrator::overlap_integral(const Orbital& o1, const Orbital& o2) const {
    return 0.0;
}

double SlaterIntegrator::electron_nucleus_integral(const Orbital& o1, const Orbital& o2, const Eigen::Vector3d& position) const {
    return 0.0;
}

double SlaterIntegrator::laplacian_integral(const Orbital& o1, const Orbital& o2) const {
    return 0.0;
}