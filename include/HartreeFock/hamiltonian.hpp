#pragma once

#include "Atom/system.hpp"
#include "Eigen/Dense"
#include "Utils/yoshimine.hpp"

class Hamiltonian {
  public:
    Hamiltonian(const System &system) { init_hamiltonian(system); }

    void init_hamiltonian(const System &system) {
        m_nuclear_repulsion = system.nuclear_energy();
        compute_one_body(system);
        compute_two_body(system);
    }

    double nuclear_repulsion() const { return m_nuclear_repulsion; }
    Eigen::MatrixXd overlap() const { return m_overlap; }
    Eigen::MatrixXd kinetic_energy() const { return m_kinetic_energy; }
    Eigen::MatrixXd electron_nuclear_energy() const {
        return m_electron_nuclear_energy;
    }
    Yoshimine<double> electron_electron_energy() const {
        return m_electron_electron_energy;
    }

  private:
    void compute_one_body(const System &system) {
        const int N = system.size();
        m_overlap.resize(N, N);
        m_kinetic_energy.resize(N, N);
        m_electron_nuclear_energy.resize(N, N);

        for (size_t j = 0; j < N; j++) {
            for (size_t i = 0; i < N; i++) {
                m_overlap(i, j) = system.overlap(i, j);
                m_overlap(j, i) = m_overlap(i, j);

                m_kinetic_energy(i, j) = system.kinetic(i, j);
                m_kinetic_energy(j, i) = m_kinetic_energy(i, j);

                m_electron_nuclear_energy(i, j) = system.electron_nucleus(i, j);
                m_electron_nuclear_energy(j, i) =
                    m_electron_nuclear_energy(i, j);
            }
        }
    }

    void compute_two_body(const System &system) {
        const int N = system.size();
        m_electron_electron_energy = Yoshimine<double>(N);

        // Yoshimine splits cases where (i > j and i < j) which are
        // equivalent,
        // we only need to compute one of them so (0 < j <= i < m) We do the
        // same for k and l. It remains cases where ij > kl and ij <= kl
        // (eg. (23|12) and (12|23) are equivalent). We only compute the
        // case where ij < kl.
        for (int i = 0; i < N; i++) {
            for (int j = 0; j <= i; j++) {
                for (int k = 0; k < N; k++) {
                    for (int l = 0; l <= k; l++) {
                        if ((i + 1) * (j + 1) > (k + 1) * (l + 1)) {
                            continue;
                        }

                        m_electron_electron_energy(i, j, k, l) =
                            system.electron_electron(i, j, k, l);
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