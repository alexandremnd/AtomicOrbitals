#pragma once

#include <vector>

#include "Atom/atom.hpp"
#include "Atom/molecule.hpp"
#include "Eigen/Dense"
#include "Integrators/electron_electron_integral.hpp"
#include "Integrators/electron_nucleus_integral.hpp"
#include "Integrators/laplacian_integral.hpp"
#include "Utils/yoshimine.hpp"
#include "Integrators/overlap_integral.hpp"

template <typename T>
class HartreeFock
{
private:
    inline void setup_system(const Atom<T>& atom) {
        m_nucleus_positions.push_back(atom.position());
        m_nucleus_charges.push_back(atom.Z());

        for (const T& orbital : atom.get_orbitals()) {
            m_orbital_basis.push_back(orbital);
        }

        m_number_of_electrons += atom.Z();
    }

    void setup_system(const std::vector<Atom<T>>& atoms) {
        m_nucleus_positions.reserve(atoms.size());
        m_nucleus_charges.reserve(atoms.size());

        for(const Atom<T>& atom : atoms) {
            setup_system(atom);
        }
    }

    void setup_overlap_matrix() {
        m_overlap = Eigen::MatrixXd::Zero(orbitals_count(), orbitals_count());

        for (size_t j = 0; j < orbitals_count(); j++) {
            for (size_t i = j; i < orbitals_count(); i++) {
                const T& orbital_i = this->orbital(i);
                const T& orbital_j = this->orbital(j);

                m_overlap(i, j) = overlap_integral(orbital_i, orbital_j);
                m_overlap(j, i) = m_overlap(i, j);
            }
        }
    }

    void setup_core_hamiltonian() {
        m_core_hamiltonian = Eigen::MatrixXd::Zero(orbitals_count(), orbitals_count());

        for (size_t j = 0; j < orbitals_count(); j++) {
            for (size_t i = j; i < orbitals_count(); i++) {
                const T& orbital_j = this->orbital(j);
                const T& orbital_i = this->orbital(i);

                m_core_hamiltonian(i, j) += -0.5 * laplacian_integral(orbital_i, orbital_j);

                for (size_t k = 0; k < atoms_count(); k++) {
                    const Eigen::Vector3d& nucleus_position = m_nucleus_positions[k];
                    const double& nucleus_charge = m_nucleus_charges[k];

                    m_core_hamiltonian(i, j) += -nucleus_charge * electron_nucleus_integral(orbital_i, orbital_j, nucleus_position);
                }

                m_core_hamiltonian(j, i) = m_core_hamiltonian(i, j);
            }
        }
    }

    void setup_electron_repulsion() {
        int N = m_orbital_basis.size();
        m_electron_repulsion = Yoshimine<double>(N);

        // Yoshimine splits cases where (i > j and i < j) which are equivalent, we only need to compute one of them so (0 < j <= i < m)
        // We do the same for k and l. It remains cases where ij > kl and ij <= kl (eg. (23|12) and (12|23) are equivalent).
        // We only compute the case where ij < kl.
        for (int i = 0; i < N; i++) {
        for (int j = 0; j <= i; j++) {
            for (int k = 0; k < N; k++) {
            for (int l = 0; l <= k; l++) {
                if ((i+1)*(j+1) > (k+1)*(l+1)) continue;

                const T& orbital_i = orbital(i);
                const T& orbital_j = orbital(j);
                const T& orbital_k = orbital(k);
                const T& orbital_l = orbital(l);
                m_electron_repulsion(i, j, k, l) = electron_electron_integral(orbital_i, orbital_j, orbital_k, orbital_l);

                std::cout << i+1 << " " << j+1 << " " << k+1 << " " << l+1 << "): " << m_electron_repulsion(i, j, k, l) << std::endl;
            }
            }
        }
        }
    }

    void diagonalize_overlap_matrix() {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(m_overlap);

        m_transformation_matrix = solver.eigenvectors() * solver.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal() * solver.eigenvectors().transpose();
    }

    void print_info() {
        std::cout << "============= Hartree-Fock Configuration =============" << std::endl;
        std::cout << "Number of atoms: " << atoms_count() << std::endl;
        std::cout << "Number of orbitals: " << orbitals_count() << std::endl;
        std::cout << "Number of electrons: " << electrons_count() << std::endl;
        std::cout << "======================================================" << std::endl;
    }

    virtual void setup_fock_matrix() = 0;
    virtual void diagonalize_fock_matrix() = 0;
    virtual void compute_density_matrix() = 0;
    virtual void compute_hf_energy() = 0;
    virtual void self_consistent_field_iteration() = 0;

protected:
    inline T& orbital(size_t i) { return m_orbital_basis[i]; }

    Eigen::MatrixXd m_overlap;
    Eigen::MatrixXd m_core_hamiltonian;
    Yoshimine<double> m_electron_repulsion;

    std::vector<T> m_orbital_basis;
    Eigen::MatrixXd m_transformation_matrix;

    std::vector<Eigen::Vector3d> m_nucleus_positions;
    std::vector<double> m_nucleus_charges;

    int m_number_of_electrons = 0;
    double m_hf_energy = 0.0;

public:
    HartreeFock(const Atom<T>& atom) {
        setup_system(atom);

        setup_overlap_matrix();
        diagonalize_overlap_matrix();

        setup_core_hamiltonian();
        setup_electron_repulsion();

        std::cout << "Overlap matrix: " << std::endl << m_overlap << std::endl;
        std::cout << "One body matrix: " << std::endl << m_core_hamiltonian << std::endl;
        std::cout << "Transformation matrix: " << std::endl << m_transformation_matrix << std::endl;
        std::cout << m_transformation_matrix.transpose() * m_overlap * m_transformation_matrix << std::endl;
        m_electron_repulsion.print_content();

        print_info();
    }

    inline int atoms_count() const { return m_nucleus_positions.size(); }
    inline int orbitals_count() const { return m_orbital_basis.size(); }
    inline int electrons_count() const { return m_number_of_electrons; }
};