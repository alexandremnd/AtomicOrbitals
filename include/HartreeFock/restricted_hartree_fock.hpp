#include "Eigen/Dense"
#include "Eigen/src/Core/Matrix.h"
#include "HartreeFock/hartree_fock.hpp"

template <typename T>
class RestrictedHartreeFock : public HartreeFock<T> {
public:
    RestrictedHartreeFock(const Atom<T>& atom) : HartreeFock<T>(atom) {
        int N = this->electrons_count();
        int L = this->orbitals_count();

        this->m_fock_matrix           = Eigen::MatrixXd::Zero(L, L);
        this->m_fock_matrix_tilde     = Eigen::MatrixXd::Zero(L, L);
        this->m_density_matrix        = Eigen::MatrixXd::Zero(L, L);
        this->m_coefficient_matrix    = Eigen::MatrixXd::Zero(L, L);

        for (int i = 0; i < 50; i++) {
            setup_fock_matrix();
            diagonalize_fock_matrix();
            compute_density_matrix();
            compute_hf_energy();

            std::cout << "HF Energy: " << this->m_hf_energy << std::endl;
        }
    }

    ~RestrictedHartreeFock() = default;

private:
    void setup_fock_matrix() override {
        int N = this->electrons_count();
        int L = this->orbitals_count();

        this->m_fock_matrix = this->m_core_hamiltonian;

        for (int p = 0; p < L; p++) {
        for (int q = 0; q < L; q++) {
            for (int r = 0; r < L; r++) {
            for (int s = 0; s < L; s++) {
                this->m_fock_matrix(q, p) += this->m_density_matrix(r, s) * this->m_electron_repulsion(p, q, r, s);
                this->m_fock_matrix(q, p) -= 0.5 * this->m_density_matrix(r, s) * this->m_electron_repulsion(p, s, r, q);
            }
            }
        }
        }

        std::cout << "Fock matrix: " << std::endl << this->m_fock_matrix << std::endl;
        std::cout << "Density matrix: " << std::endl << this->m_density_matrix << std::endl;
    }

    void diagonalize_fock_matrix() override {
        int N = this->electrons_count();
        int L = this->orbitals_count();

        const Eigen::MatrixXd &F = m_fock_matrix;
        const Eigen::MatrixXd &V = this->m_transformation_matrix;
        Eigen::MatrixXd &F_tilde = this->m_fock_matrix_tilde;
        Eigen::MatrixXd &C = this->m_coefficient_matrix;


        std::cout << "\n";
        std::cout << "Fock matrix: " << std::endl << F << std::endl;
        std::cout << "Transformation matrix: " << std::endl << V << std::endl;

        F_tilde = V * F * V;

        std::cout <<  "Transformed Fock matrix: " << std::endl << F_tilde << std::endl;

        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(F_tilde);

        std::cout << "Eigenvalues: " << std::endl << es.eigenvalues() << std::endl;
        std::cout << "Eigenvectors: " << std::endl << es.eigenvectors() << std::endl;

        C = V * es.eigenvectors().block(0, 0, L, N/2);
        std::cout << "Coefficient matrix: " << std::endl << C << std::endl;
    }

    void compute_density_matrix() override {
        this->m_density_matrix = 0.3 * 2 * this->m_coefficient_matrix * this->m_coefficient_matrix.transpose() + 0.7 * this->m_density_matrix;
    }

    void self_consistent_field_iteration() override {
        this->setup_fock_matrix();
        this->diagonalize_fock_matrix();
        this->compute_density_matrix();
        this->compute_hf_energy();
    }

    void compute_hf_energy() override {
        int N = this->electrons_count();
        int L = this->orbitals_count();

        this->m_hf_energy = 0.0;

        for (int p = 0; p < L; p++) {
        for (int q = 0; q < L; q++) {
            this->m_hf_energy += 0.5 * this->m_density_matrix(p, q) * (this->m_core_hamiltonian(p, q) + this->m_fock_matrix(p, q));
        }
        }
    }

    Eigen::MatrixXd m_fock_matrix;
    Eigen::MatrixXd m_fock_matrix_tilde;
    Eigen::MatrixXd m_density_matrix;
    Eigen::MatrixXd m_coefficient_matrix;
};