#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE

#include "include/Eigen/Dense"
#include <iostream>

#include "include/HartreeFock/hartree_fock.hpp"
#include "include/BasisSet/gaussian_primitive.hpp"

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

void copyright_warranty() {
    std::cout << "======================================\n";
    std::cout << "AtomicOrbitals Copyright (C) 2024 by Ewan Bataille, Jéremy Atané and Alexandre Menard\n";
    std::cout << "This program comes with ABSOLUTELY NO WARRANTY;\n";
    std::cout << "This is free software, and you are welcome to\n";
    std::cout << "redistribute it under certain conditions;\n";
    std::cout << "======================================\n\n\n";
}

int main() {
    copyright_warranty();
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