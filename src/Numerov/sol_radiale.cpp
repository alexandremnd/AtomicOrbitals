// Import de librairies
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>
#include "Numerov/numerov.hpp"
#include "Numerov/sol_radiale.hpp"

void radial_solution(int Z, int l, int n) {
    const double E = -Z * Z / (2. * n * n); // Energie de l'électron

    double dr = 0.001, ri = 200. / Z,
           rf = 0.001; // Le rayon initiale ne doit pas être trop grand pour
                       // éviter une divergence en 0. Le rayon initiale doit
                       // prendre en compte que si Z grand, la zone de
                       // l'électron est plus proche du noyau donc les
                       // conditions initiales sont à mettre plus proche aussi
    double Nr = fabs(rf - ri) / dr; // Nombre de points de calcul
    std::vector<double> Y0(2), R0(2);
    Y0[0] = 0.0001;
    Y0[1] = 0.0001;
    R0[0] = ri;
    R0[1] = ri - dr;

    auto Ec = [E, Z, l](double r) { return Energie(r, E, Z, l); };

    std::vector<double> U = Numerov(Y0, Ec, f_nulle, R0, ri, rf, dr);

    std::string Nom_fichier = "valeur_radiale.txt";
    std::ofstream fichier(Nom_fichier);
    fichier << "E=" << E << ";Z=" << Z << ";l=" << l << ";n=" << n
            << ";dr=" << dr << std::endl;
    fichier << "Résultat de la méthode de Numérov avec les valeurs de la "
               "fonction, les valeurs en r"
            << std::endl;
    fichier << "y(r);r" << std::endl;

    double yj, rj;
    for (int j = 0; j < Nr; j++) {
        yj = U[j];
        rj = ri - dr * j;
        fichier << std::to_string(yj) << ";" << std::to_string(rj) << std::endl;
    }
    fichier.close();

    std::cout << "dr=" << dr << ", Nr=" << Nr << std::endl;
    std::cout << "E=" << E << " Hartree" << std::endl;
    std::cout << "Z=" << Z << ", l=" << l << ", n=" << n << std::endl;
}