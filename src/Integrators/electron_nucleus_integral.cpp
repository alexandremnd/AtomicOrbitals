#include "Integrators/electron_nucleus_integral.hpp"
#include <stdexcept>


double electron_nucleus_integral(const GaussianPrimitive& orbital1, const GaussianPrimitive& orbital2, const Eigen::Vector3d& nucleus_position) {
    throw std::logic_error("Not implemented yet");
}


double electron_nucleus_integral(const SlaterPrimitive& orbital1, const SlaterPrimitive& orbital2, const Eigen::Vector3d& nucleus_position) {
    throw std::logic_error("Not implemented yet");
}