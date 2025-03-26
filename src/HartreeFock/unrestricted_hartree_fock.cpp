#include "HartreeFock/unrestricted_hartree_fock.hpp"
#include <cstddef>
#include <iostream>

UnrestrictedHartreeFock::UnrestrictedHartreeFock(const System &system,
                                                 uint num_electrons,
                                                 uint num_alpha_electrons)
    : HartreeFock(system, num_electrons),
      m_num_alpha_electrons(num_alpha_electrons),
      m_num_beta_electrons(num_electrons - num_alpha_electrons) {

    const size_t L = m_H.size();
    const size_t N = m_no_electrons;

    // Initialize alpha and beta specific matrices
    m_alpha_density_matrix = Eigen::MatrixXd::Zero(L, L);
    m_beta_density_matrix = Eigen::MatrixXd::Zero(L, L);
    m_alpha_fock_matrix = Eigen::MatrixXd::Zero(L, L);
    m_beta_fock_matrix = Eigen::MatrixXd::Zero(L, L);
    m_alpha_fock_matrix_tilde = Eigen::MatrixXd::Zero(L, L);
    m_beta_fock_matrix_tilde = Eigen::MatrixXd::Zero(L, L);
    m_alpha_coefficients = Eigen::MatrixXd::Zero(L, m_num_alpha_electrons);
    m_beta_coefficients = Eigen::MatrixXd::Zero(L, m_num_beta_electrons);

    std::cout << "UnrestrictedHartreeFock initialized with "
              << m_num_alpha_electrons << " alpha electrons and "
              << m_num_beta_electrons << " beta electrons." << std::endl;
}

void UnrestrictedHartreeFock::setup_fock_matrix() {
    int N = m_no_electrons;
    int L = m_H.size();

    // Initialize Fock matrices with the core Hamiltonian
    m_alpha_fock_matrix = m_H.T() + m_H.V();
    m_beta_fock_matrix = m_H.T() + m_H.V();

    for (size_t p = 0; p < L; p++) {
        for (size_t q = 0; q < L; q++) {
            for (size_t r = 0; r < L; r++) {
                for (size_t s = 0; s < L; s++) {
                    double pqrs = m_H.ee(p, q, r, s);
                    double prsq = m_H.ee(p, r, q, s);
                    m_alpha_fock_matrix(p, q) +=
                        (pqrs - prsq) * m_alpha_density_matrix(r, s) +
                        pqrs * m_beta_density_matrix(r, s);
                    m_beta_fock_matrix(p, q) +=
                        (pqrs - prsq) * m_beta_density_matrix(r, s) +
                        pqrs * m_alpha_density_matrix(r, s);
                }
            }
        }
    }
}

void UnrestrictedHartreeFock::diagonalize_fock_matrix() {
    int N = m_no_electrons;
    int L = m_H.size();

    const Eigen::MatrixXd &V = m_transformation_matrix;
    Eigen::MatrixXd &F_alpha_tilde = m_alpha_fock_matrix_tilde;
    Eigen::MatrixXd &C_alpha = m_alpha_coefficients;
    Eigen::MatrixXd &F_beta_tilde = m_beta_fock_matrix_tilde;
    Eigen::MatrixXd &C_beta = m_beta_coefficients;

    F_alpha_tilde = V.transpose() * m_alpha_fock_matrix * V;

    F_beta_tilde = V.transpose() * m_beta_fock_matrix * V;

    // Solve eigenvalue problems for alpha and beta electrons
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> alpha_solver(F_alpha_tilde);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> beta_solver(F_beta_tilde);

    // Transform eigenvectors back to original basis
    m_alpha_coefficients =
        V * alpha_solver.eigenvectors().block(0, 0, L, m_num_alpha_electrons);
    m_beta_coefficients =
        V * beta_solver.eigenvectors().block(0, 0, L, m_num_beta_electrons);

    // Normalize coefficient matrices
    normalizeCoefficientMatrix();
}

void UnrestrictedHartreeFock::normalizeCoefficientMatrix() {
    const size_t L = m_H.size();
    const size_t N = m_no_electrons;

    // Normalize alpha coefficient matrix
    for (size_t k = 0; k < m_num_alpha_electrons; k++) {
        double normalizationFactor = 0;
        for (size_t p = 0; p < L; p++) {
            for (size_t q = 0; q < L; q++) {
                normalizationFactor += m_alpha_coefficients(p, k) *
                                       m_alpha_coefficients(q, k) *
                                       m_H.S()(p, q);
            }
        }
        normalizationFactor = sqrt(normalizationFactor);
        m_alpha_coefficients.col(k) =
            m_alpha_coefficients.col(k) / normalizationFactor;
    }

    // Normalize beta coefficient matrix
    for (size_t k = 0; k < m_num_beta_electrons; k++) {
        double normalizationFactor = 0;
        for (size_t p = 0; p < L; p++) {
            for (size_t q = 0; q < L; q++) {
                normalizationFactor += m_beta_coefficients(p, k) *
                                       m_beta_coefficients(q, k) *
                                       m_H.S()(p, q);
            }
        }
        normalizationFactor = sqrt(normalizationFactor);
        m_beta_coefficients.col(k) =
            m_beta_coefficients.col(k) / normalizationFactor;
    }
}

void UnrestrictedHartreeFock::compute_density_matrix() {
    const size_t L = m_H.size();
    const size_t N = m_no_electrons;

    // Reset density matrices
    m_alpha_density_matrix = m_smoothing_factor * m_alpha_coefficients *
                                 m_alpha_coefficients.transpose() +
                             (1 - m_smoothing_factor) * m_alpha_density_matrix;
    m_beta_density_matrix = m_smoothing_factor * m_beta_coefficients *
                                m_beta_coefficients.transpose() +
                            (1 - m_smoothing_factor) * m_beta_density_matrix;
}

void UnrestrictedHartreeFock::compute_hf_energy() {
    double electronic_energy = 0.0;

    electronic_energy +=
        (m_H.T() + m_H.V())
            .cwiseProduct(m_alpha_density_matrix + m_beta_density_matrix)
            .sum();

    electronic_energy +=
        m_alpha_fock_matrix.cwiseProduct(m_alpha_density_matrix).sum();
    electronic_energy +=
        m_beta_fock_matrix.cwiseProduct(m_beta_density_matrix).sum();

    electronic_energy /= 2.0;

    m_hf_energy = electronic_energy + m_H.nucleus_repulsion();
}

void UnrestrictedHartreeFock::self_consistent_field_iteration(
    size_t iteration) {
    setup_fock_matrix();
    diagonalize_fock_matrix();
    compute_density_matrix();

    compute_hf_energy();
}