// Import de librairies
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "Numerov/numerov.hpp"
#include "Numerov/sol_radiale.hpp"

void radial_solution(int Z, int l, int n) {
    const double E = -Z * Z / (2. * n * n); // Energie de l'électron

    // Le rayon initiale ne doit pas être trop grand pour
    // éviter une divergence en 0. Le rayon initiale doit
    // prendre en compte que si Z grand, la zone de
    // l'électron est plus proche du noyau donc les
    // conditions initiales sont à mettre plus proche aussi

    double dr = 0.001;
    double ri = 200. / Z;
    double rf = 0.001;

    double Nr = fabs(rf - ri) / dr; // Nombre de points de calcul

    Eigen::VectorXd Y0(2), R0(2);

    // Conditions initiales
    Y0[0] = 0.0001;
    Y0[1] = 0.0001;
    R0[0] = ri;
    R0[1] = ri - dr;

    auto Ec = [E, Z, l](double r) { return Numerov::energy(r, E, Z, l); };

    Eigen::VectorXd U =
        Numerov::numerov(Y0, Ec, Numerov::f_null, R0, ri, rf, dr);

    std::string output_file = "data/valeur_radiale.out";
    std::ofstream fichier(output_file);

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

    // Affichage des paramètres
    std::cout << "dr=" << dr << ", Nr=" << Nr << std::endl;
    std::cout << "E=" << E << " Hartree" << std::endl;
    std::cout << "Z=" << Z << ", l=" << l << ", n=" << n << std::endl;
}