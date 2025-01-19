#include "HartreeFock/hartree_fock.hpp"

template <typename T>
class RestrictedHartreeFock : public HartreeFock<T> {
public:
    RestrictedHartreeFock(const Atom<T>& atom) : HartreeFock<T>(atom) {
        m_fock_matrix = Eigen::MatrixXd::Zero(this->m_orbital_basis.size(), this->m_orbital_basis.size());
        m_density_matrix = Eigen::MatrixXd::Zero(this->m_orbital_basis.size(), this->m_orbital_basis.size());
    }

    RestrictedHartreeFock(const Molecule<T>& molecule) : HartreeFock<T>(molecule) {
        m_fock_matrix = Eigen::MatrixXd::Zero(this->m_orbital_basis.size(), this->m_orbital_basis.size());
        m_density_matrix = Eigen::MatrixXd::Zero(this->m_orbital_basis.size(), this->m_orbital_basis.size());
    }

    ~RestrictedHartreeFock() = default;

private:
    void setup_fock_matrix() override {
        const int N = this->m_orbital_basis.size();

        for (int j = 0; j < N; j++) {
        for (int i = j; i < N; i++) {
            this->m_fock_matrix(i, j) += this->m_one_body_matrix(i, j);

            for (int k = 0; k < N; k++) {
            for (int l = 0; l < N; l++) {
                const T& orbital_i = this->orbital(i);
                const T& orbital_j = this->orbital(j);
                const T& orbital_k = this->orbital(k);
                const T& orbital_l = this->orbital(l);

                this->m_fock_matrix(i, j) += this->m_density_matrix(k, l) * this->m_two_body_matrix(i, j, k, l);
                this->m_fock_matrix(i, j) -= 0.5 * this->m_density_matrix(k, l) * this->m_two_body_matrix(i, l, k, j);
            }
            }

            this->m_fock_matrix(j, i) = this->m_fock_matrix(i, j);
        }
        }
    }

    void diagonalize_fock_matrix() override {

    }

    void compute_density_matrix() override {

    }

    void compute_hf_energy() override {

    }

private:
    Eigen::MatrixXd m_fock_matrix;
    Eigen::MatrixXd m_density_matrix;
    double m_hf_energy;
};