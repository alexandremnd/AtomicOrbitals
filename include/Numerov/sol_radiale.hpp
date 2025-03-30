#pragma once

/**
 * @brief Computes the radial solution for an atomic orbital
 *
 * Calculates the radial wavefunction for an atomic orbital using the Numerov
 * method to solve the radial Schrödinger equation for a hydrogen-like atom.
 *
 * @param Z The atomic number of the element
 * @param l The angular momentum quantum number
 * @param n The principal quantum number
 *
 * @note This function may store the results internally or output them to a file
 * @see Numerov algorithm for differential equations
 */
void radial_solution(int Z, int l, int n);