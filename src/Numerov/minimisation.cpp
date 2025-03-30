// Import de librairies
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>

#include "Numerov/numerov.hpp"
#include "Numerov/minimisation.hpp"

std::vector<double> Multiplication(std::vector<double> U,
                                   std::vector<double> V) {
    int Taille = U.size();
    std::vector<double> W(Taille);
    for (int i = 0; i < Taille; i++) {
        W[i] = U[i] * V[i];
    }
    return W;
}

std::vector<double> Multiplication_scalaire(std::vector<double> U, double x) {
    int Taille = U.size();
    std::vector<double> W(Taille);
    for (int i = 0; i < Taille; i++) {
        W[i] = U[i] * x;
    }
    return W;
}

double maximum(std::vector<double> L) {
    double maxi = L[0];
    for (int i = 0; i < L.size(); i++) {
        if (L[i] > maxi) {
            maxi = L[i];
        }
    }
    return maxi;
}

std::vector<double> derivation(std::vector<double> &X, std::vector<double> &Y) {
    int Taille = X.size();
    double dx, dy;
    std::vector<double> dY(Taille - 1);
    for (int i = 0; i < Taille - 1; i++) {
        dx = X[i + 1] - X[i];
        dy = Y[i + 1] - Y[i];
        dY[i] = dy / dx;
    }
    return dY;
}

/// Fonctions utiles

std::vector<double>
Liste_Energie(double Ei, double Ef,
              double dE_i) { // Fonctionne uniquement pour des énergies
                             // négatives tq |Ei| > |Ef|
    std::vector<double> L_E;
    double E = Ei, dE = dE_i;

    while (E <= Ef) {
        if (fabs(E / Ef) < 10 && fabs(E / dE) < 100) {
            dE = dE / 10.;
        }
        L_E.push_back(E);
        E = E + dE;
    }

    L_E.shrink_to_fit(); // Fixe la taille du vecteur une fois initialisé
    return L_E;
}

std::vector<double> f_U0(double dr, double ri, double rf, int Nr,
                         std::vector<double> &L_E, int N_E, int Z, int l) {
    std::vector<double> Y0(2), r0(2), U(Nr), L_U0(N_E), densite;

    Y0[0] = 0.0001;
    Y0[1] = 0.0001;
    r0[0] = ri;
    r0[1] = ri - dr;

    double E, max_dense;

    for (int i = 0; i < N_E; i++) { // On parcourt les énergies
        E = L_E[i];
        auto Ec = [E, Z, l](double r) { return Energie(r, E, Z, l); };
        U = Numerov(Y0, Ec, f_nulle, r0, ri, rf, dr);

        densite = Multiplication(U, Multiplication_scalaire(U, dr));
        max_dense = maximum(densite);
        L_U0[i] = (densite[Nr - 1]) / max_dense;
    }
    return L_U0;
}

void minimisation(int Z, int l) {
    double
        dr = 0.001,
        ri = 200. /
             Z, // Le rayon initiale ne doit pas être trop grand pour éviter une
                // divergence en 0. Le rayon initiale doit prendre en compte que
                // si Z grand, la zone de l'électron est plus proche du noyau
                // donc les conditions initiales sont à mettre plus proche aussi
        rf = 0.01; // Le rayon finale, soit le plus proche de 0, doit être
                   // adapté à la valeur de l
    double Nr = fabs(rf - ri) / dr;
    double dE = 0.0002, Ei = -0.7,
           Ef = -0.05; // Les énergies doivent être négatives et |Ei| > |Ef|
    std::vector<double> L_E = Liste_Energie(Ei, Ef, dE);
    int N_E = L_E.size();

    std::cout << "dr=" << dr << ", Nr=" << Nr << std::endl;
    std::cout << "dE=" << dE << ", N_E=" << N_E << std::endl;
    std::vector<double> L_U0 = f_U0(dr, ri, rf, Nr, L_E, N_E, Z, l);
    std::vector<double> L_U0_norm =
        Multiplication_scalaire(L_U0, 1. / maximum(L_U0));

    std::string Nom_fichier = "valeur_rfixe.txt"; // Nom du fichier de sortie
    std::ofstream fichier(Nom_fichier);
    fichier << "Z=" << Z << ";l=" << l << ";dr=" << dr << std::endl;
    fichier << "U(0,E);E" << std::endl;

    double L_U0j, Ej;
    for (int j = 0; j < N_E; j++) {
        L_U0j = L_U0_norm[j];
        Ej = L_E[j];
        fichier << std::to_string(L_U0j) << ";" << std::to_string(Ej)
                << std::endl;
    }
    fichier.close();
}