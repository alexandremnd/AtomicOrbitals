#pragma once

#include "Orbitals/gaussian_primitive.hpp"
#include "Eigen/src/Core/Matrix.h"
#include "Maths/tensor4D.hpp"

#include <cmath>

// Compute the coefficient for the decomposition of the product of two
// gaussian on the base of hermite gaussian

bool ConditionRecurrence(int n, int t, int u, int v);

Tensor4D<double> HermiteIntegral(const GaussianPrimitive &orbital1,
                                 const GaussianPrimitive &orbital2, double p,
                                 Eigen::Vector3d position);