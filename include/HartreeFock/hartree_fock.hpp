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
    }

    void setup_system(const std::vector<Atom<T>>& atoms) {
        m_nucleus_positions.reserve(atoms.size());
        m_nucleus_charges.reserve(atoms.size());

        for(const Atom<T>& atom : atoms) {
            setup_system(atom);
        }
    }

    void setup_overlap_matrix() {
        m_overlap_matrix = Eigen::MatrixXd::Zero(orbitals_count(), orbitals_count());

        for (size_t j = 0; j < orbitals_count(); j++) {
            for (size_t i = j; i < orbitals_count(); i++) {
                const T& orbital_i = this->orbital(i);
                const T& orbital_j = this->orbital(j);

                m_overlap_matrix(i, j) = overlap_integral(orbital_i, orbital_j);
                m_overlap_matrix(j, i) = m_overlap_matrix(i, j);
            }
        }
    }

    void setup_one_body_integrals() {
        m_one_body_matrix = Eigen::MatrixXd::Zero(orbitals_count(), orbitals_count());

        for (size_t j = 0; j < orbitals_count(); j++) {
            for (size_t i = j; i < orbitals_count(); i++) {
                const T& orbital_j = this->orbital(j);
                const T& orbital_i = this->orbital(i);

                m_one_body_matrix(i, j) += -0.5 * laplacian_integral(orbital_i, orbital_j);

                for (size_t k = 0; k < atoms_count(); k++) {
                    m_one_body_matrix(i, j) -= m_nucleus_charges[k] * electron_nucleus_integral(orbital_i, orbital_j, m_nucleus_positions[k]);
                }

                m_one_body_matrix(j, i) = m_one_body_matrix(i, j);
            }
        }
    }

    void setup_two_body_integrals() {
        int N = m_orbital_basis.size();
        m_two_body_matrix = Yoshimine<double>((N*(N+1)*(N*N+N+2))/8);

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
                m_two_body_matrix(i, j, k, l) = electron_electron_integral(orbital_i, orbital_j, orbital_k, orbital_l);
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

protected:
    inline T& orbital(size_t i) { return m_orbital_basis[i]; }

    Eigen::MatrixXd m_overlap_matrix;
    Eigen::MatrixXd m_one_body_matrix;
    Yoshimine<double> m_two_body_matrix;

    std::vector<T> m_orbital_basis;
    Eigen::MatrixXd m_transformation_matrix;

    std::vector<Eigen::Vector3d> m_nucleus_positions;
    std::vector<double> m_nucleus_charges;

public:
    HartreeFock(const Atom<T>& atom) {
        setup_system(atom);

        setup_overlap_matrix();
        diagonalize_overlap_matrix();

        setup_one_body_integrals();
        setup_two_body_integrals();
    }

    HartreeFock(const Molecule<T>& molecule) {
        setup_system(molecule.atoms());

        setup_overlap_matrix();
        diagonalize_overlap_matrix();

        setup_one_body_integrals();
        setup_two_body_integrals();
    }

    inline int atoms_count() const { return m_nucleus_positions.size(); }
    inline int orbitals_count() const { return m_orbital_basis.size(); }
};