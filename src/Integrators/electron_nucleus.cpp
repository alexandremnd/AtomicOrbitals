#include "Integrators/electron_nucleus.hpp"
#include "BasisSet/gaussian_primitive.hpp"
#include "BasisSet/slater_primitive.hpp"
#include "Integrators/overlap_primitive.hpp"
#include <stdexcept>


double electron_nucleus_integral(const GaussianPrimitive& orbital1, const GaussianPrimitive& orbital2, const Eigen::Vector3d& nucleus_position) {
    throw std::logic_error("Not implemented yet");
}

double electron_nucleus_integral(const SlaterPrimitive& orbital1, const SlaterPrimitive& orbital2, const Eigen::Vector3d& nucleus_position) {
    assert(nucleus_position.isZero()); // Slater HF computation are not implemented for multi-atoms systems

    return overlap_integral(orbital1, orbital2, -1);
}