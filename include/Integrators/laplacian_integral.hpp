#pragma once

#include "include/BasisSet/gaussian_primitive.hpp"
#include "include/BasisSet/slater_primitive.hpp"

template <typename T>
double laplacian_integral(const T& orbital1, const T& orbital2);

template <>
double laplacian_integral(const GaussianPrimitive& orbital1,
                        const GaussianPrimitive& orbital2);

template <>
double laplacian_integral(const SlaterPrimitive& orbital1,
                        const SlaterPrimitive& orbital2);
