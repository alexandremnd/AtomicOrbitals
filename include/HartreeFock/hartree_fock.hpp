#pragma once

#include "Atom/system.hpp"
#include "Eigen/Dense"
#include "HartreeFock/hamiltonian.hpp"

class HartreeFock {
  private:
    void diagonalize_overlap_matrix();

    virtual void setup_fock_matrix() = 0;
    virtual void diagonalize_fock_matrix() = 0;
    virtual void compute_density_matrix() = 0;
    virtual void compute_hf_energy() = 0;
    virtual void self_consistent_field_iteration() = 0;

  protected:
    Eigen::MatrixXd m_transformation_matrix;
    Hamiltonian m_H;

    double m_hf_energy = 1000.0;
    int m_no_electrons;

  public:
    HartreeFock(const System &system, uint no_electrons) {
        set_system(system);
        m_no_electrons = no_electrons;

        print_info();
    }

    void set_system(const System &system);
    void run();

    void print_info();
    void print_iteration_info(int n);
    void print_result(int total_iterations);
};