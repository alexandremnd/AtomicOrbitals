#pragma once

#include "BasisSet/gaussian_primitive.hpp"
#include "BasisSet/gaussian_contracted.hpp"
#include "BasisSet/slater_primitive.hpp"

double laplacian_integral(const GaussianPrimitive& orbital1,
                        const GaussianPrimitive& orbital2);

double laplacian_integral(const ContractedGaussian& orbital1,
                        const ContractedGaussian& orbital2);

double laplacian_integral(const SlaterPrimitive& orbital1,
                        const SlaterPrimitive& orbital2);
