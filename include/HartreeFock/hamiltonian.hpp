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
        return 2 * m_electron_electron_energy(i, j, k, l) -
               m_electron_electron_energy(i, l, k, j);
    }

    size_t size() const { return m_overlap.rows(); }

    void print() const {
        std::cout << "Overlap matrix:\n" << m_overlap << std::endl;
        for (size_t i = 0; i < size(); i++) {
            for (size_t j = 0; j < size(); j++) {
                std::cout << i << " " << j << " " << m_overlap(i, j)
                          << std::endl;
            }
        }

        std::cout << "Kinetic energy matrix:\n";
        for (size_t i = 0; i < size(); i++) {
            for (size_t j = 0; j < size(); j++) {
                std::cout << i << " " << j << " " << m_kinetic_energy(i, j)
                          << std::endl;
            }
        }

        std::cout << "Potential matrix \n";
        for (size_t i = 0; i < size(); i++) {
            for (size_t j = 0; j < size(); j++) {
                std::cout << i << " " << j << " "
                          << m_electron_nuclear_energy(i, j) << std::endl;
            }
        }

        std::cout << "Electron-electron energy tensor\n";
        for (size_t i = 0; i < size(); i++) {
            for (size_t j = 0; j < size(); j++) {
                for (size_t k = 0; k < size(); k++) {
                    for (size_t l = 0; l < size(); l++) {
                        std::cout << i << " " << j << " " << k << " " << l
                                  << " "
                                  << m_electron_electron_energy(i, j, k, l)
                                  << std::endl;
                    }
                }
            }
        }
    }

  private:
    void compute_one_body(const System &system);
    void compute_two_body(const System &system);

    double m_nuclear_repulsion;
    Eigen::MatrixXd m_overlap;
    Eigen::MatrixXd m_kinetic_energy;
    Eigen::MatrixXd m_electron_nuclear_energy;
    Yoshimine<double> m_electron_electron_energy;
};