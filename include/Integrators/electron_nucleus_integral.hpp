#pragma once

#include "include/BasisSet/gaussian_primitive.hpp"
#include "include/BasisSet/slater_primitive.hpp"

#include "include/Eigen/Dense"

template <typename T>
double electron_nucleus_integral(const T& orbital1, const T& orbital2, const Eigen::Vector3d& nucleus_position);

template <>
double electron_nucleus_integral(const GaussianPrimitive& orbital1,
                                const GaussianPrimitive& orbital2,
                                const Eigen::Vector3d& nucleus_position);

template <>
double electron_nucleus_integral(const SlaterPrimitive& orbital1,
                                const SlaterPrimitive& orbital2,
                                const Eigen::Vector3d& nucleus_position);
