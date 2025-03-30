#pragma once

#include <iostream>

typedef unsigned long size_t;

/**
 * @brief Any system of particles passed to Hartree-Fock must implement this
 * class.
 *
 */
class System {
  public:
    /**
     * @brief Overlap integral between two basis functions i and j
     *
     * @param i i-th basis function
     * @param j j-th basis function
     * @return double <i|j>
     */
    virtual double overlap(size_t i, size_t j) const = 0;

    /**
     * @brief Kinetic energy integral between two basis functions i and j
     *
     * @param i i-th basis function
     * @param j j-th basis function
     * @return double -0.5 <i|∇²|j>
     */
    virtual double kinetic(size_t i, size_t j) const = 0;

    /**
     * @brief Electron-nucleus integral between two basis functions i and j
     *
     * @param i i-th basis function
     * @param j j-th basis function
     * @return double -Z <i| 1/|r-rA| |j>
     */
    virtual double electron_nucleus(size_t i, size_t j) const = 0;

    /**
     * @brief Electron-electron integral between four basis functions i, j, k, l
     *
     * @param i i-th basis function
     * @param j j-th basis function sharing same coordinates as i
     * @param k k-th basis function
     * @param l l-th basis function sharing same coordinates as k
     * @return double <ik|w|jl> := (ij|kl)
     */
    virtual double electron_electron(size_t i, size_t j, size_t k,
                                     size_t l) const = 0;

    /**
     * @brief Repulsion energy between all nuclei
     *
     * @return double Sum of all repulsion energies between nuclei
     */
    virtual double nucleus_repulsion() const = 0;

    virtual size_t size() const = 0;
};