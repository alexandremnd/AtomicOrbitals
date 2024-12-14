#pragma once

#include "BasisSet/slater_primitive.hpp"
#include "BasisSet/gaussian_primitive.hpp"
#include "BasisSet/gaussian_contracted.hpp"

double overlap_integral(const GaussianPrimitive& orbital1, const GaussianPrimitive& orbital2);

double overlap_integral(const ContractedGaussian& orbital1, const ContractedGaussian& orbital2);

double overlap_integral(const SlaterPrimitive& orbital1, const SlaterPrimitive& orbital2, const int n_offset = 0);