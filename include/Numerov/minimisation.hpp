#pragma once

#include "Eigen/Dense"

/**
 * @brief Performs energy minimization for a given atomic system using the
 * Numerov method
 *
 * @param Z Atomic number of the system
 * @param l Angular momentum quantum number
 */
void minimisation(int Z, int l);

/**
 * @brief Calculates the numerical derivative of a function represented by (X,Y)
 * data points
 *
 * @param X Vector of x-coordinates
 * @param Y Vector of y-coordinates (function values at X)
 * @return std::vector<double> Vector containing the numerical derivatives of Y
 * with respect to X
 */
Eigen::VectorXd derivative(Eigen::VectorXd &X, Eigen::VectorXd &Y);

/**
 * @brief Generates a list of energy values from Ei to Ef with specified initial
 * step (and decreasing step size logarithmically)
 * Works only for negative energies such that |Ei| > |Ef|.
 *
 * @param Ei Initial energy value
 * @param Ef Final energy value
 * @param dE_i Initial energy step size
 * @return std::vector<double> List of energy values from Ei to Ef
 */
Eigen::VectorXd energy_list(double Ei, double Ef, double dE_i);

/**
 * @brief Computes initial wavefunctions using Numerov's method for a given set
 * of parameters
 *
 * @param dr Radial step size
 * @param ri Initial radius value
 * @param rf Final radius value
 * @param Nr Number of radial points
 * @param L_E List of energy values to evaluate
 * @param N_E Number of energy values
 * @param Z Atomic number
 * @param l Angular momentum quantum number
 * @return std::vector<double> Resulting wavefunction values
 */
Eigen::VectorXd f_U0(double dr, double ri, double rf, int Nr,
                     Eigen::VectorXd &L_E, int Z, int l);