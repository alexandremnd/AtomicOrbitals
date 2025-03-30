#pragma once

#include <vector>
#include <functional>

/**
 * @brief Implements the Numerov method to solve second-order differential
 * equations of the form y''(x) = -g(x)y(x) + s(x).
 *
 * The Numerov method is a numerical technique commonly used in quantum
 * mechanics and other fields to solve differential equations to a high degree
 * of accuracy.
 *
 * @param Y0 Reference to vector containing initial values of the solution.
 * @param f_g Function representing g(x) in the differential equation.
 * @param f_s Function representing s(x) in the differential equation.
 * @param X0 Reference to vector containing initial grid points.
 * @param xi Starting point of the integration range.
 * @param xf Ending point of the integration range.
 * @param dx Step size for integration.
 *
 * @return Vector containing the solution values at the specified grid points.
 */
std::vector<double> Numerov(std::vector<double> &Y0,
                            std::function<double(double)> f_g,
                            double (&f_s)(double), std::vector<double> &X0,
                            double xi, double xf, double dx);

double V_eff(double r, int Z, int l);

double Energie(double r, double E, int Z, int l);

double f_nulle(double x);