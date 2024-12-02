#pragma once

#include "include/BasisSet/basis_set.hpp"

template <typename T>
double overlap_integral(const T& orbital1, const T& orbital2);

template <>
double overlap_integral(const GaussianPrimitive& orbital1, const GaussianPrimitive& orbital2);

template <>
double overlap_integral(const SlaterPrimitive& orbital1, const SlaterPrimitive& orbital2);