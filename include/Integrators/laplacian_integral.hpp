#pragma once

#include "include/BasisSet/gaussian_primitive.hpp"
#include "include/BasisSet/slater_primitive.hpp"

double laplacian_integral(const GaussianPrimitive& orbital1,
                        const GaussianPrimitive& orbital2);

double laplacian_integral(const SlaterPrimitive& orbital1,
                        const SlaterPrimitive& orbital2);
