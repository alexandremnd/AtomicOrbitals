#pragma once

#include "Orbitals/gaussian_primitive.hpp"
#include "Eigen/src/Core/Matrix.h"
#include "Maths/tensor4D.hpp"

#include <cmath>

bool ConditionRecurrence(int n, int t, int u, int v);

/**
 * @brief Computes the coefficient for the decomposition of the product of two
 * gaussians on the base of hermite gaussians integral
 *
 * @param orbital1 First orbital to be decomposed
 * @param orbital2 Second orbital to be decomposed
 * @param p
 * @param position Position of the nucleus
 * @return Tensor4D<double>
 */
Tensor4D<double> HermiteIntegral(const GaussianPrimitive &orbital1,
                                 const GaussianPrimitive &orbital2, double p,
                                 Eigen::Vector3d position);