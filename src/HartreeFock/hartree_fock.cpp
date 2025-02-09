#include <iomanip>
#include <iostream>

#include "Eigen/Dense"
#include "HartreeFock/hartree_fock.hpp"

void HartreeFock::diagonalize_overlap_matrix() {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(m_H.S());

    m_transformation_matrix =
        solver.eigenvectors() *
        solver.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal() *
        solver.eigenvectors().transpose();
}

void HartreeFock::set_system(const System &system) {
    m_H = Hamiltonian(system);
    diagonalize_overlap_matrix();
}

void HartreeFock::set_smoothing_factor(double smoothing_factor) {
    if (smoothing_factor < 0 || smoothing_factor > 1) {
        throw std::invalid_argument("Smoothing factor must be between 0 and 1");
    }
    m_smoothing_factor = smoothing_factor;
}

void HartreeFock::run(double convergence_threshold, int max_iterations,
                      bool silent) {
    size_t n = 0;
    double old_energy;

    do {
        old_energy = m_hf_energy;
        self_consistent_field_iteration(n);
        if (!silent) {
            print_iteration_info(n);
        }
        n++;
    } while (std::abs(m_hf_energy - old_energy) > convergence_threshold &&
             n < max_iterations);

    print_result(n);
}

void HartreeFock::print_info() {
    std::cout << "============= Hartree-Fock Configuration ============="
              << std::endl;
    std::cout << "Number of electrons: " << m_no_electrons << std::endl;
    std::cout << "Number of basis functions: " << m_H.size() << std::endl;
    std::cout << "======================================================"
              << std::endl;
}

void HartreeFock::print_iteration_info(int n) {
    std::cout << std::setprecision(8);
    std::cout << "Iteration n°" << n << ": E(SCF) = " << m_hf_energy
              << " (Hartree)\n";
}

void HartreeFock::print_result(int total_iterations) {
    std::cout << "============= Hartree-Fock result ===================="
              << std::endl;
    std::cout << "SCF Iteration required: " << total_iterations << "\n";
    std::cout << "Final energy: " << m_hf_energy << " (Hartree)\n";
    std::cout << "Final energy: " << m_hf_energy * 27.2113 << " (eV)\n";
}