#pragma once

#include "Eigen/Dense"

#include <vector>
#include <functional>

namespace Numerov {
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
Eigen::VectorXd numerov(Eigen::VectorXd &Y0, std::function<double(double)> f_g,
                        double (&f_s)(double), Eigen::VectorXd &X0, double xi,
                        double xf, double dx);

/**
 * @brief Effective potential for the radial Schrödinger equation.
 *
 * @param r Radial position
 * @param Z Nucleus charge
 * @param l Angular momentum quantum number
 * @return double Effective potential value
 */
double v_eff(double r, int Z, int l);

/**
 * @brief Energy function for the radial Schrödinger equation.
 *
 * @param r Radial position
 * @param E Energy
 * @param Z Nucleus charge
 * @param l Agular momentum quantum number
 * @return double
 */
double energy(double r, double E, int Z, int l);

/**
 * @brief Function that returns 0.
 *
 * @param x Input value
 * @return double Always returns 0.
 */
double f_null(double x);

} // namespace Numerov