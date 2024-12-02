#include "include/Integrators/laplacian_integral.hpp"
#include <stdexcept>

template <>
double laplacian_integral(const GaussianPrimitive& orbital1, const GaussianPrimitive& orbital2) {
    throw std::logic_error("Not implemented yet");
}

template <>
double laplacian_integral(const SlaterPrimitive& orbital1, const SlaterPrimitive& orbital2) {
    throw std::logic_error("Not implemented yet");
}