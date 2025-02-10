#pragma once

#include "BasisSet/gaussian_primitive.hpp"
#include "Maths/tensor3D.hpp"

#include <cmath>

// Compute the coefficient for the decomposition of the product of two gaussian
// on the base of hermite gaussian

bool ConditionIndex(int i, int j, int t);

std::vector<Tensor3D<double>>
HermiteCoefficient(const GaussianPrimitive &orbital1,
                   const GaussianPrimitive &orbital2);