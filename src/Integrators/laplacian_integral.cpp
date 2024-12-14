#include "Integrators/laplacian_integral.hpp"
#include "BasisSet/slater_primitive.hpp"
#include "Integrators/overlap_integral.hpp"
#include <stdexcept>

double laplacian_integral(const GaussianPrimitive& orbital1, const GaussianPrimitive& orbital2) {
    throw std::logic_error("Not implemented yet");
}

double laplacian_integral(const ContractedGaussian& orbital1, const ContractedGaussian& orbital2) {
    throw std::logic_error("Not implemented yet");
}

double laplacian_integral(const SlaterPrimitive& orbital1, const SlaterPrimitive& orbital2) {
    double matrix_element = orbital2.alpha() * orbital2.alpha() * overlap_integral(orbital1, orbital2);
    matrix_element -= 2 * orbital2.alpha() * orbital2.n() * overlap_integral(orbital1, orbital2, -1);
    matrix_element += (orbital2.n() * (orbital2.n() - 1) - orbital2.l() * (orbital2.l() + 1)) * overlap_integral(orbital1, orbital2, -2);

    return matrix_element;
}