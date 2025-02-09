#include "HartreeFock/restricted_hartree_fock.hpp"
#include "Eigen/Dense"

RestrictedHartreeFock::RestrictedHartreeFock(const System &system,
                                             uint no_electrons)
    : HartreeFock(system, no_electrons) {
    int N = m_no_electrons;
    int L = m_H.size();

    m_fock_matrix = Eigen::MatrixXd::Zero(L, L);
    m_fock_matrix_tilde = Eigen::MatrixXd::Zero(L, L);
    m_density_matrix = Eigen::MatrixXd::Zero(L, L);
    m_coefficient_matrix = Eigen::MatrixXd::Zero(L, N / 2);
}

void RestrictedHartreeFock::setup_fock_matrix() {
    int N = m_no_electrons;
    int L = m_H.size();

    m_fock_matrix = m_H.T() + m_H.V();

    for (int p = 0; p < L; p++) {
        for (int q = 0; q < L; q++) {
            for (int r = 0; r < L; r++) {
                for (int s = 0; s < L; s++) {
                    m_fock_matrix(q, p) +=
                        0.5 * m_density_matrix(r, s) * m_H.as(p, q, r, s);
                }
            }
        }
    }
}

void RestrictedHartreeFock::normalizeCoefficientMatrix() {
    int N = m_no_electrons;
    int L = m_H.size();

    for (int k = 0; k < m_no_electrons / 2; k++) {
        double normalizationFactor = 0;
        for (int p = 0; p < L; p++) {
            for (int q = 0; q < L; q++) {
                normalizationFactor += m_coefficient_matrix(p, k) *
                                       m_coefficient_matrix(q, k) *
                                       m_H.S()(p, q);
            }
        }
        normalizationFactor = sqrt(normalizationFactor);
        m_coefficient_matrix.col(k) =
            m_coefficient_matrix.col(k) / normalizationFactor;
    }
}

void RestrictedHartreeFock::diagonalize_fock_matrix() {
    int N = m_no_electrons;
    int L = m_H.size();

    const Eigen::MatrixXd &F = m_fock_matrix;
    const Eigen::MatrixXd &V = m_transformation_matrix;
    Eigen::MatrixXd &F_tilde = m_fock_matrix_tilde;
    Eigen::MatrixXd &C = m_coefficient_matrix;

    F_tilde = V * F * V;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(F_tilde);

    C = V * es.eigenvectors().block(0, 0, L, N / 2);
    normalizeCoefficientMatrix();
}

void RestrictedHartreeFock::compute_density_matrix() {
    m_density_matrix = m_smoothing_factor * 2 * m_coefficient_matrix *
                           m_coefficient_matrix.transpose() +
                       (1 - m_smoothing_factor) * m_density_matrix;
}

void RestrictedHartreeFock::self_consistent_field_iteration() {
    setup_fock_matrix();
    diagonalize_fock_matrix();
    compute_density_matrix();
    compute_hf_energy();
    std::cout << m_density_matrix << std::endl;
}

void RestrictedHartreeFock::compute_hf_energy() {
    int N = m_no_electrons;
    int L = m_H.size();

    m_hf_energy = 0;

    for (int p = 0; p < L; p++)
        for (int q = 0; q < L; q++) {
            m_hf_energy += m_density_matrix(p, q) * m_H.core(p, q);

            for (int r = 0; r < L; r++)
                for (int s = 0; s < L; s++) {
                    m_hf_energy += m_H.as(p, q, r, s) * m_density_matrix(p, q) *
                                   m_density_matrix(s, r) * 0.25;
                }
        }
    m_hf_energy += m_H.nuclear_repulsion();
}
