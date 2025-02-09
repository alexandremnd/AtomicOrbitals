#pragma once

#include "Eigen/Dense"
#include "HartreeFock/hartree_fock.hpp"

class RestrictedHartreeFock : public HartreeFock {
  public:
    RestrictedHartreeFock(const System &system, uint no_electrons);
    ~RestrictedHartreeFock() = default;

  private:
    void setup_fock_matrix() override;
    void diagonalize_fock_matrix() override;
    void compute_density_matrix() override;
    void self_consistent_field_iteration() override;
    void compute_hf_energy() override;

    void normalizeCoefficientMatrix();

    Eigen::MatrixXd m_fock_matrix;
    Eigen::MatrixXd m_fock_matrix_tilde;
    Eigen::MatrixXd m_density_matrix;
    Eigen::MatrixXd m_coefficient_matrix;
};