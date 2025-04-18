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

#pragma omp parallel for
    for (size_t p = 0; p < L; p++) {
        for (size_t q = 0; q < L; q++) {
            for (size_t r = 0; r < L; r++) {
                for (size_t s = 0; s < L; s++) {
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

    for (size_t k = 0; k < m_no_electrons / 2; k++) {
        double normalizationFactor = 0;
        for (size_t p = 0; p < L; p++) {
            for (size_t q = 0; q < L; q++) {
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
    Eigen::VectorXd &E = m_orbital_energies;

    F_tilde = V * F * V;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(F_tilde);

    E = es.eigenvalues().head(N / 2);
    C = V * es.eigenvectors().block(0, 0, L, N / 2);
    normalizeCoefficientMatrix();
}

void RestrictedHartreeFock::compute_density_matrix() {
    m_density_matrix = m_smoothing_factor * 2 * m_coefficient_matrix *
                           m_coefficient_matrix.transpose() +
                       (1 - m_smoothing_factor) * m_density_matrix;
}

void RestrictedHartreeFock::self_consistent_field_iteration(size_t iteration) {
    setup_fock_matrix();
    diagonalize_fock_matrix();
    compute_density_matrix();

    if (m_diis_size > 0) {
        diis(iteration);
    }

    compute_hf_energy();
}

void RestrictedHartreeFock::compute_hf_energy() {
    int N = m_no_electrons;
    int L = m_H.size();

    m_hf_energy =
        0.5 *
        m_density_matrix.cwiseProduct(m_fock_matrix + m_H.T() + m_H.V()).sum();

    m_hf_energy += m_H.nucleus_repulsion();
}

void RestrictedHartreeFock::reset_diis_subspace() {
    m_fock_history.clear();
    m_density_history.clear();
    m_error_history =
        Eigen::MatrixXd::Zero(m_H.size() * m_H.size(), m_diis_size);
    m_fock_history.resize(m_diis_size);
    m_density_history.resize(m_diis_size);
}

void RestrictedHartreeFock::diis(size_t iteration) {
    Eigen::MatrixXd error_vector =
        m_fock_matrix * m_density_matrix - m_density_matrix * m_fock_matrix;
    error_vector.resize(error_vector.size(), 1);

    // We update the DIIS (direct inversion in iterative subspace) error history
    size_t index = iteration % m_diis_size;
    m_error_history.col(index) = error_vector;
    m_fock_history[index] = m_fock_matrix;
    m_density_history[index] = m_density_matrix;

    if (iteration > m_diis_size) {
        diis_compute();
        return;
    }

    if (iteration % (2 * m_diis_size) == 0) {
        reset_diis_subspace();
    }
}

void RestrictedHartreeFock::diis_compute() {
    Eigen::MatrixXd gram_matrix =
        Eigen::MatrixXd::Zero(m_diis_size + 1, m_diis_size + 1);

    gram_matrix.block(0, 0, m_diis_size, m_diis_size) =
        m_error_history.transpose() * m_error_history;

    gram_matrix.row(m_diis_size).setConstant(1);
    gram_matrix.col(m_diis_size).setConstant(-1);
    gram_matrix(m_diis_size, m_diis_size) = 0;

    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(m_diis_size + 1);
    rhs(m_diis_size) = 1;

    Eigen::VectorXd ci = gram_matrix.colPivHouseholderQr().solve(rhs);

    m_fock_matrix = Eigen::MatrixXd::Zero(m_H.size(), m_H.size());
    m_density_matrix = Eigen::MatrixXd::Zero(m_H.size(), m_H.size());

    for (size_t i = 0; i < m_diis_size; i++) {
        m_fock_matrix += ci(i) * m_fock_history[i];
        m_density_matrix += ci(i) * m_density_history[i];
    }
}