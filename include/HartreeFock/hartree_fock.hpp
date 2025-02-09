#pragma once

#include <iomanip>
#include <iostream>

#include "Eigen/Dense"
#include "concepts.hpp"

#include "HartreeFock/hartree_fock.interface.hpp" // IWYU pragma: export

template <DerivedFromOrbital OrbitalType>
void HartreeFock<OrbitalType>::diagonalize_overlap_matrix() {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(m_overlap);

    m_transformation_matrix =
        solver.eigenvectors() *
        solver.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal() *
        solver.eigenvectors().transpose();
}

template <DerivedFromOrbital OrbitalType>
void HartreeFock<OrbitalType>::print_info() {
    std::cout << "============= Hartree-Fock Configuration ============="
              << std::endl;
    std::cout << "Number of electrons: " << electrons_count() << std::endl;
    std::cout << "======================================================"
              << std::endl;
}

template <DerivedFromOrbital OrbitalType>
void HartreeFock<OrbitalType>::print_iteration_info(int n) {
    std::cout << std::setprecision(8);
    std::cout << "Iteration n°" << n << ": E(SCF) = " << m_hf_energy
              << " (Hartree)\n";
}

template <DerivedFromOrbital OrbitalType>
void HartreeFock<OrbitalType>::print_result(int total_iterations) {
    std::cout << "======================================================"
              << std::endl;
    std::cout << "SCF Iteration required: " << total_iterations << "\n";
    std::cout << "Final energy: " << m_hf_energy << " (Hartree)\n";
    std::cout << "Final energy: " << m_hf_energy * 27.2113 << " (eV)\n";
}