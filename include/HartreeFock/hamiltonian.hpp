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
    void compute_one_body(const System &system);
    void compute_two_body(const System &system);

    double m_nuclear_repulsion;
    Eigen::MatrixXd m_overlap;
    Eigen::MatrixXd m_kinetic_energy;
    Eigen::MatrixXd m_electron_nuclear_energy;
    Yoshimine<double> m_electron_electron_energy;
};