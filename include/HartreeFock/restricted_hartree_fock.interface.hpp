#pragma once

#include "Eigen/Dense"
#include "Atom/atom.hpp"
#include "HartreeFock/hartree_fock.hpp"
#include "concepts.hpp"


template <DerivedFromOrbital OrbitalType>
class RestrictedHartreeFock : public HartreeFock<OrbitalType> {
public:
    RestrictedHartreeFock(const Atom<OrbitalType>& atom);
    ~RestrictedHartreeFock() = default;

private:
    void setup_fock_matrix() override;
    void diagonalize_fock_matrix() override;
    void compute_density_matrix() override;
    void self_consistent_field_iteration() override;
    void compute_hf_energy() override;

    Eigen::MatrixXd m_fock_matrix;
    Eigen::MatrixXd m_fock_matrix_tilde;
    Eigen::MatrixXd m_density_matrix;
    Eigen::MatrixXd m_coefficient_matrix;
};

DECLARE_EXTERN_TEMPLATE(RestrictedHartreeFock)