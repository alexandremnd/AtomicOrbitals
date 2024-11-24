#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE

#include "include/Eigen/Dense"
#include <iostream>

Eigen::MatrixXd potential_matrix(int N, double r_min, double dr, int L, int Z) {
    Eigen::MatrixXd V(N, N);
    for (int i = 0; i < N; i++) {
        double r = r_min + i * dr;
            V(i, i) = L * (L+1) / (r*r) - 2 * Z / r;
    }
    return V;
}

Eigen::MatrixXd kinetic_matrix(int N, double dr) {
    Eigen::MatrixXd T(N, N);
    for (int i = 0; i < N; i++) {
        T(i, i) = -2/(dr * dr);
        if (i < N-1) {
            T(i, i+1) = 1/(dr*dr);
            T(i+1, i) = 1/(dr*dr);
        }
    }
    return T;
}

int main() {
    double r_min = 0.00001;
    double r_max = 50.0;
    int N = 1000;
    int L = 0;
    int Z = 1;
    double dr = (r_max - r_min) / N;

    Eigen::MatrixXd hamiltonian = -kinetic_matrix(N, dr) + potential_matrix(N, r_min, dr, L, Z);

    // eigenvalues of selfadjoint m
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(hamiltonian);

    for (size_t i = 0; i < 10; i++)
    {
        std::cout << "The " << i << "th eigenvalue is: " << es.eigenvalues()[i]/2 << std::endl;
    }


    return 0;
}