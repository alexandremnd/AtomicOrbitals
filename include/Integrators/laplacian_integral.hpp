#pragma once

#include "BasisSet/slater_primitive.hpp"
#include "BasisSet/gaussian_primitive.hpp"
#include "BasisSet/contracted_orbital.hpp"

double laplacian_integral(const GaussianPrimitive& orbital1,
                        const GaussianPrimitive& orbital2);

double laplacian_integral(const ContractedGaussian& orbital1,
                        const ContractedGaussian& orbital2);

double laplacian_integral(const SlaterPrimitive& orbital1,
                        const SlaterPrimitive& orbital2);
