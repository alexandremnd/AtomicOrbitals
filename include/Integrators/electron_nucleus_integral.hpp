#pragma once

#include "include/BasisSet/gaussian_primitive.hpp"
#include "include/BasisSet/gaussian_contracted.hpp"
#include "include/BasisSet/slater_primitive.hpp"

#include "include/Eigen/Dense"

double electron_nucleus_integral(const GaussianPrimitive& orbital1,
                                const GaussianPrimitive& orbital2,
                                const Eigen::Vector3d& nucleus_position);

double electron_nucleus_integral(const ContractedGaussian& orbital1,
                                const ContractedGaussian& orbital2,
                                const Eigen::Vector3d& nucleus_position);

double electron_nucleus_integral(const SlaterPrimitive& orbital1,
                                const SlaterPrimitive& orbital2,
                                const Eigen::Vector3d& nucleus_position);
