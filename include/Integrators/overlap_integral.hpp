#pragma once

#include "include/BasisSet/basis_set.hpp"

double overlap_integral(const GaussianPrimitive& orbital1, const GaussianPrimitive& orbital2);

double overlap_integral(const SlaterPrimitive& orbital1, const SlaterPrimitive& orbital2);