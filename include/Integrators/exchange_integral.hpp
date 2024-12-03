#pragma once

#include "include/BasisSet/gaussian_primitive.hpp"
#include "include/BasisSet/slater_primitive.hpp"


double exchange_integral(const GaussianPrimitive& orbital1,
                                const GaussianPrimitive& orbital2,
                                const GaussianPrimitive& orbital3,
                                const GaussianPrimitive& orbital4);

double exchange_integral(const SlaterPrimitive& orbital1,
                                const SlaterPrimitive& orbital2,
                                const SlaterPrimitive& orbital3,
                                const SlaterPrimitive& orbital4);
