#pragma once

#include "Atom/system.hpp"
#include "Eigen/Dense"
#include "Utils/yoshimine.hpp"

class Hamiltonian {
  public:
    Hamiltonian() = default;
    Hamiltonian(const System &system) { compute_hamiltonian(system); }

    void compute_hamiltonian(const System &system) {
        m_nuclear_repulsion = system.nuclear_energy();
        compute_one_body(system);
        compute_two_body(system);
    }

    /**
     * @brief Compute the nuclear repulsion energy between all nuclei in the
     * system.
     */
    double nuclear_repulsion() const { return m_nuclear_repulsion; }

    /**
     * @brief Returns the overlap matrix.
     *
     * @return const Eigen::MatrixXd& Overlap matrix.
     */
    const Eigen::MatrixXd &S() const { return m_overlap; }

    /**
     * @brief Returns the kinetic energy matrix.
     *
     * @return const Eigen::MatrixXd& Kinetic energy matrix.
     */
    const Eigen::MatrixXd &T() const { return m_kinetic_energy; }

    /**
     * @brief Returns the electron-nucleus attraction matrix.
     *
     * @return const Eigen::MatrixXd& Electron-nucleus attraction matrix.
     */
    const Eigen::MatrixXd &V() const { return m_electron_nuclear_energy; }

    /**
     * @brief Returns the electron-electron repulsion tensor.
     *
     * @return const Yoshimine<double>& Electron-electron repulsion tensor.
     */
    const Yoshimine<double> &EE() const { return m_electron_electron_energy; }

    /**
     * @brief Returns the core Hamiltonian matrix.
     *
     * @param i i-th orbital
     * @param j j-th orbital
     * @return double One body term.
     */
    double core(size_t i, size_t j) const {
        return m_kinetic_energy(i, j) + m_electron_nuclear_energy(i, j);
    }

    /**
     * @brief Returns the antisymmetric two-electron integral.
     *
     * @param i i-th orbital
     * @param j j-th orbital
     * @param k k-th orbital
     * @param l l-th orbital
     * @return double
     */
    double as(size_t i, size_t j, size_t k, size_t l) const {
        return m_electron_electron_energy(i, j, k, l) -
               0.5 * m_electron_electron_energy(i, l, k, j);
    }

    size_t size() const { return m_overlap.rows(); }

  private:
    void compute_one_body(const System &system);
    void compute_two_body(const System &system);

    double m_nuclear_repulsion;
    Eigen::MatrixXd m_overlap;
    Eigen::MatrixXd m_kinetic_energy;
    Eigen::MatrixXd m_electron_nuclear_energy;
    Yoshimine<double> m_electron_electron_energy;
};