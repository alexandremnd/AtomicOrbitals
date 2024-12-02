#include "include/Integrators/overlap_integral.hpp"
#include <stdexcept>

template <>
double overlap_integral(const GaussianPrimitive& orbital1, const GaussianPrimitive& orbital2) {
    throw std::logic_error("Not implemented");
}

template <>
double overlap_integral(const SlaterPrimitive& orbital1, const SlaterPrimitive& orbital2) {
    throw std::logic_error("Not implemented");
}
