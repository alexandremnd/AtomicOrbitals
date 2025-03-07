#pragma once

#include "Eigen/Dense"
#include "Eigen/src/Core/Matrix.h"
#include "HartreeFock/hartree_fock.hpp"
#include <sys/types.h>

/**
 * @brief Hartree-Fock method implementation in the restricted case. Basis
 * orbitals are supposed doubly occupied.
 *
 */
class RestrictedHartreeFock : public HartreeFock {
  public:
    RestrictedHartreeFock(const System &system, uint no_electrons);
    ~RestrictedHartreeFock() = default;

    void set_diis_size(uint size) {
        m_diis_size = size;
        reset_diis_subspace();
    }

    void reset_diis_subspace() {
        m_fock_history.clear();
        m_density_history.clear();
        m_error_history =
            Eigen::MatrixXd::Zero(m_H.size() * m_H.size(), m_diis_size);
        m_fock_history.resize(m_diis_size);
        m_density_history.resize(m_diis_size);
    }

  private:
    void setup_fock_matrix() override;
    void diagonalize_fock_matrix() override;
    void compute_density_matrix() override;
    void self_consistent_field_iteration(size_t iteration) override;
    void compute_hf_energy() override;

    void normalizeCoefficientMatrix();

    /**
     * @brief Stores the error vector in the DIIS error history
     *
     * @param iteration Iteration count of the SCF loop. Used to determine the
     * index of the error vector in the error history
     */
    void diis(size_t iteration);

    /**
     * @brief Computes the DIIS coefficients and updates the Fock and density
     * matrices accordingly.
     * @see https://en.wikipedia.org/wiki/DIIS
     */
    void diis_compute();

    Eigen::MatrixXd m_fock_matrix;
    Eigen::MatrixXd m_fock_matrix_tilde;
    Eigen::MatrixXd m_density_matrix;
    Eigen::MatrixXd m_coefficient_matrix;

    uint m_diis_size = 0;
    std::vector<Eigen::MatrixXd> m_fock_history;
    std::vector<Eigen::MatrixXd> m_density_history;
    Eigen::MatrixXd m_error_history;
};