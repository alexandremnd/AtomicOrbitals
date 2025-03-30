#pragma once

#include "HartreeFock/hartree_fock.hpp"
#include <Eigen/Dense>

class UnrestrictedHartreeFock : public HartreeFock {
  public:
    UnrestrictedHartreeFock(const System &system, uint num_electrons,
                            uint num_alpha_electrons);

    const Eigen::VectorXd &orbital_alpha_energies() {
        return m_alpha_orbital_energies;
    }
    const Eigen::MatrixXd &coefficient_alpha_matrix() {
        return m_alpha_coefficients;
    }
    const Eigen::MatrixXd &density_alpha_matrix() {
        return m_alpha_density_matrix;
    }

    const Eigen::VectorXd &orbital_beta_energies() {
        return m_beta_orbital_energies;
    }

    const Eigen::MatrixXd &coefficient_beta_matrix() {
        return m_beta_coefficients;
    }

    const Eigen::MatrixXd &density_beta_matrix() {
        return m_beta_density_matrix;
    }

  private:
    void setup_fock_matrix() override;
    void diagonalize_fock_matrix() override;
    void compute_density_matrix() override;
    void self_consistent_field_iteration(size_t iteration) override;
    void compute_hf_energy() override;

    void normalizeCoefficientMatrix();

    int m_num_alpha_electrons;
    int m_num_beta_electrons;
    Eigen::MatrixXd m_alpha_density_matrix;
    Eigen::MatrixXd m_beta_density_matrix;
    Eigen::MatrixXd m_alpha_fock_matrix;
    Eigen::MatrixXd m_beta_fock_matrix;
    Eigen::MatrixXd m_alpha_fock_matrix_tilde;
    Eigen::MatrixXd m_beta_fock_matrix_tilde;
    Eigen::MatrixXd m_alpha_coefficients;
    Eigen::MatrixXd m_beta_coefficients;
    Eigen::VectorXd m_alpha_orbital_energies;
    Eigen::VectorXd m_beta_orbital_energies;
};
