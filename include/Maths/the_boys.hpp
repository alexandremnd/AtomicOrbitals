#pragma once

#include <boost/math/special_functions/gamma.hpp>
#include <cmath>

/**
 * @brief Computes the Boys function of order n evaluated at x
 *
 * @param n Order of the function
 * @param x Point of evaluation
 * @return double Evaluated function Fn(x)
 */
double TheBoys_Function(int n, double x);

/**
 * @brief Computes recursively the values of the boys function from order O to
 * order n.
 *
 * @param n Maximum order to be computed
 * @param x Point of evaluation
 * @return std::vector<double> Values of Fm(x) for 0<=m<=n : [F0,F1,...,Fn]
 */
std::vector<double> TheBoys_Recurrence(int n, double x);