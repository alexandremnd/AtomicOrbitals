// Import de librairies
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>

#include "Eigen/Dense"
#include "Numerov/numerov.hpp"
#include "Numerov/minimisation.hpp"

Eigen::VectorXd derivative(Eigen::VectorXd &X, Eigen::VectorXd &Y) {
    int N = X.size();
    Eigen::VectorXd dx, dy;

    Eigen::VectorXd dY(N - 1);

    // dx = X[1:] - X[:-1]
    // dy = Y[1:] - Y[:-1]
    dx = X.block(1, 0, N - 1, 1) - X.block(0, 0, N - 1, 1);
    dy = Y.block(1, 0, N - 1, 1) - Y.block(0, 0, N - 1, 1);

    // dY = dy / dx
    dY = dy.cwiseQuotient(dx);

    return dY;
}

Eigen::VectorXd energy_list(double Ei, double Ef, double dE_i) {
    std::vector<double> L_E;
    double E = Ei, dE = dE_i;

    while (E <= Ef) {
        if (fabs(E / Ef) < 10 && fabs(E / dE) < 100) {
            dE = dE / 10.;
        }
        L_E.push_back(E);
        E = E + dE;
    }

    L_E.shrink_to_fit();

    Eigen::VectorXd out(L_E.size());
    for (size_t i = 0; i < L_E.size(); i++) {
        out(i) = L_E[i];
    }

    return out;
}

Eigen::VectorXd f_U0(double dr, double ri, double rf, int Nr,
                     Eigen::VectorXd &L_E, int Z, int l) {
    int N_E = L_E.size();
    Eigen::VectorXd Y0(2), r0(2), U(Nr), L_U0(N_E), densite;

    // Initialisation des conditions initiales
    Y0[0] = 0.0001;
    Y0[1] = 0.0001;
    r0[0] = ri;
    r0[1] = ri - dr;

    double E, max_dense;

    for (int i = 0; i < N_E; i++) {
        E = L_E[i];

        auto Ec = [E, Z, l](double r) { return Numerov::energy(r, E, Z, l); };
        U = Numerov::numerov(Y0, Ec, Numerov::f_null, r0, ri, rf, dr);

        // On calcule la densité normalisée du dernier point (Nr - 1) ie en 0
        densite = U.cwiseProduct(U * dr); // U * U * dr
        max_dense = densite.maxCoeff();
        L_U0[i] = (densite[Nr - 1]) / max_dense;
    }

    return L_U0;
}

void minimisation(int Z, int l) {
    // Le rayon initial (ri) ne doit pas être trop grand pour éviter une
    // divergence en 0. Le rayon initiale doit prendre en compte que
    // si Z grand, la zone de l'électron est plus proche du noyau
    // donc les conditions initiales sont à mettre plus proche aussi

    // Le rayon final doit être le plus proche de 0 et doit être
    // adapté à la valeur de l

    double dr = 0.001;
    double ri = 200. / Z;
    double rf = 0.01;

    double Nr = fabs(rf - ri) / dr;
    double dE = 0.0002;
    double Ei = -0.7;
    double Ef = -0.05; // Les énergies doivent être négatives et |Ei| > |Ef|

    Eigen::VectorXd L_E = energy_list(Ei, Ef, dE);
    int N_E = L_E.size();

    // Affichage des paramètres
    std::cout << "dr=" << dr << ", Nr=" << Nr << std::endl;
    std::cout << "dE=" << dE << ", N_E=" << N_E << std::endl;

    Eigen::VectorXd L_U0 = f_U0(dr, ri, rf, Nr, L_E, Z, l);
    Eigen::VectorXd L_U0_norm = L_U0 / L_U0.maxCoeff();

    std::string output_file = "data/valeur_rfixe.out";
    std::ofstream fichier(output_file);

    fichier << "Z=" << Z << ";l=" << l << ";dr=" << dr << std::endl;
    fichier << "U(0,E);E" << std::endl;

    for (int j = 0; j < N_E; j++) {
        fichier << std::to_string(L_U0_norm[j]) << ";" << std::to_string(L_E[j])
                << std::endl;
    }

    fichier.close();
}