#pragma once

#include "Eigen/Dense"
#include "Utils/yoshimine.hpp"
#include "concepts.hpp"

template <DerivedFromOrbital OrbitalType>
class HartreeFock {
private:
    void diagonalize_overlap_matrix();

    void print_info();
    void print_iteration_info(int n);
    void print_result(int total_iterations);

    virtual void setup_fock_matrix() = 0;
    virtual void diagonalize_fock_matrix() = 0;
    virtual void compute_density_matrix() = 0;
    virtual void compute_hf_energy() = 0;
    virtual void self_consistent_field_iteration() = 0;

protected:
    Eigen::MatrixXd m_overlap;
    Eigen::MatrixXd m_core_hamiltonian;
    Yoshimine<double> m_electron_repulsion;

    Eigen::MatrixXd m_transformation_matrix;

    int m_number_of_electrons = 0;
    int m_number_of_orbitals = 0;
    double m_hf_energy = 0.0;

public:
    HartreeFock() {
        diagonalize_overlap_matrix();

        print_info();
    }

    inline int electrons_count() const { return m_number_of_electrons; }
    inline int orbitals_count() const { return m_number_of_orbitals; }
};