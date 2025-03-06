#pragma once

#include "Orbitals/gaussian_primitive.hpp"
#include "Maths/tensor3D.hpp"

#include <cmath>

bool ConditionIndex(int i, int j, int t);

/**
 * @brief Computesthe coefficient for the decomposition of the product of two
 * gaussians on the basis of hermite gaussians
 *
 * @param orbital1 First orbital to be decomposed
 * @param orbital2 Second orbital to be decomposed
 * @return std::vector<Tensor3D<double>> Coefficients for the decomposition for
 * each axis (x, y, z)
 */
std::vector<Tensor3D<double>>
HermiteCoefficient(const GaussianPrimitive &orbital1,
                   const GaussianPrimitive &orbital2);