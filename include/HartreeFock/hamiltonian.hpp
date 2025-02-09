#pragma once

#include "Atom/molecule.hpp"
#include "BasisSet/slater_primitive.hpp"
#include "Eigen/Dense"
#include "Utils/yoshimine.hpp"
#include "concepts.hpp"

template <DerivedFromOrbital OrbitalType> class Hamiltonian {
  public:
    Hamiltonian(const Atom<OrbitalType> &atom) {
        compute_one_body(atom.get_orbitals(), {Eigen::Vector3d::Zero()},
                         {atom.Z()});
        compute_two_body(atom.get_orbitals());
    }

    Hamiltonian(const Molecule<OrbitalType> &molecule) {
        compute_one_body(molecule.get_orbitals(), molecule.get_positions(),
                         molecule.get_charges());
        compute_two_body(molecule.orbitals());

        for (size_t i = 0; i < molecule.size(); i++) {
            for (size_t j = i + 1; j < molecule.size(); j++) {
                m_nuclear_repulsion += molecule.get_charges[i] *
                                       molecule.get_charges[j] *
                                       molecule.distance(i, j);
            }
        }
    }

  private:
    void compute_one_body(std::vector<OrbitalType> &orbitals,
                          std::vector<Eigen::Vector3d> &positions,
                          std::vector<int> &charges) {

        const int N = orbitals.size();
        m_overlap.resize(N, N);
        m_kinetic_energy.resize(N, N);
        m_electron_nuclear_energy.resize(N, N);

        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < N; j++) {
                m_overlap(j, i) = overlap_integral(orbitals[j], orbitals[j]);
                m_overlap(i, j) = m_overlap(j, i);

                m_kinetic_energy(j, i) =
                    -0.5 * laplacian_integral(orbitals[j], orbitals[i]);
                m_kinetic_energy(i, j) = m_kinetic_energy(j, i);

                for (size_t k = 0; k < positions.size(); k++) {
                    m_electron_nuclear_energy(j, i) +=
                        -charges[k] * electron_nuclear_integral(orbitals[j],
                                                                orbitals[i],
                                                                positions[k]);
                }
            }
        }
    }

    void compute_two_body(std::vector<OrbitalType> &orbitals) {
        const int N = orbitals.size();
        m_electron_electron_energy = Yoshimine<double>(N);

        // Yoshimine splits cases where (i > j and i < j) which are equivalent,
        // we only need to compute one of them so (0 < j <= i < m) We do the
        // same for k and l. It remains cases where ij > kl and ij <= kl (eg.
        // (23|12) and (12|23) are equivalent). We only compute the case where
        // ij < kl.
        for (int i = 0; i < N; i++) {
            for (int j = 0; j <= i; j++) {
                for (int k = 0; k < N; k++) {
                    for (int l = 0; l <= k; l++) {
                        if ((i + 1) * (j + 1) > (k + 1) * (l + 1)) {
                            continue;
                        }

                        const OrbitalType &orbital_i = orbitals[i];
                        const OrbitalType &orbital_j = orbitals[j];
                        const OrbitalType &orbital_k = orbitals[k];
                        const OrbitalType &orbital_l = orbitals[l];
                        m_electron_electron_energy(i, j, k, l) =
                            electron_electron_integral(orbital_i, orbital_j,
                                                       orbital_k, orbital_l);
                    }
                }
            }
        }
    }

    double m_nuclear_repulsion;
    Eigen::MatrixXd m_overlap;
    Eigen::MatrixXd m_kinetic_energy;
    Eigen::MatrixXd m_electron_nuclear_energy;
    Yoshimine<double> m_electron_electron_energy;
};