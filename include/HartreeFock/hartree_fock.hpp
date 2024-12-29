#pragma once

#include <iostream>
#include <vector>

#include "Atom/atom.hpp"
#include "Atom/molecule.hpp"
#include "Eigen/Dense"
#include "Integrators/electron_electron_integral.hpp"
#include "Utils/yoshimine.hpp"
#include "Integrators/overlap_integral.hpp"

template <typename T>
class HartreeFock
{
private:
    void setup_basis(const Atom<T>& atom) {
        for (const T& orbital : atom.get_orbitals()) {
            m_orbital_basis.push_back(orbital);
        }
    }

    void setup_basis(const std::vector<Atom<T>>& atoms) {
        for (const Atom<T>& atom : atoms) {
            setup_basis(atom);
        }
    }

    void setup_overlap_matrix() {
        m_overlap_matrix = Eigen::MatrixXd::Zero(m_orbital_basis.size(), m_orbital_basis.size());

        for (size_t j = 0; j < m_orbital_basis.size(); j++) {
            for (size_t i = j; i < m_orbital_basis.size(); i++) {
                m_overlap_matrix(i, j) = overlap_integral(m_orbital_basis[i], m_orbital_basis[j]);
                m_overlap_matrix(j, i) = m_overlap_matrix(i, j);
            }
        }
    }

    void setup_two_body_integrals() {
        int m = m_orbital_basis.size();
        m_yoshimine = Yoshimine<double>((m*(m+1)*(m*m+m+2))/8);

        // Yoshimine split cases where (i > j and i < j) which are equivalent, we only need to compute one of them so (0 < j <= i < m)
        // We do the same for k and l. It remains cases where ij > kl and ij <= kl (eg. (23|12) and (12|23) are equivalent).
        // We only compute the case where ij < kl.
        for (int i = 0; i < m; i++) {
        for (int j = 0; j <= i; j++) {
            for (int k = 0; k < m; k++) {
            for (int l = 0; l <= k; l++) {
                if ((i+1)*(j+1) > (k+1)*(l+1)) continue;

                const T& orbital_i = m_orbital_basis[i];
                const T& orbital_j = m_orbital_basis[j];
                const T& orbital_k = m_orbital_basis[k];
                const T& orbital_l = m_orbital_basis[l];
                m_yoshimine(i, j, k, l) = electron_electron_integral(orbital_i, orbital_j, orbital_k, orbital_l);
            }
            }
        }
        }
    }

    void diagonalize_overlap_matrix() {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(m_overlap_matrix);

        m_transformation_matrix = solver.eigenvectors() * solver.eigenvalues().cwiseSqrt().cwiseInverse().asDiagonal();
    }

    virtual void setup_fock_matrix() = 0;
    virtual void diagonalize_fock_matrix() = 0;
    virtual void compute_density_matrix() = 0;
    virtual void compute_hf_energy() = 0;

    Eigen::MatrixXd m_transformation_matrix;
    Eigen::MatrixXd m_overlap_matrix;
    std::vector<T> m_orbital_basis;
    Yoshimine<double> m_yoshimine;

public:
    HartreeFock(const Atom<T>& atom) {
        setup_basis(atom);
        setup_overlap_matrix();
        diagonalize_overlap_matrix();
        setup_two_body_integrals();
    }

    HartreeFock(const Molecule<T>& molecule) {
        setup_basis(molecule.atoms());
        setup_overlap_matrix();
        diagonalize_overlap_matrix();
        setup_two_body_integrals();
    }
};