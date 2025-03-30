// Import de librairies
#include <iostream>
#include <vector>
#include <functional>
#include <math.h>

#include "Numerov/numerov.hpp"

double V_eff(double r, int Z, int l) {
    return -Z / r + l * (l + 1) / (2. * r * r);
}

double Energie(double r, double E, int Z, int l) {
    return 2. * (E - V_eff(r, Z, l));
}

double f_nulle(double x) { return 0.; }

std::vector<double> Numerov(std::vector<double> &Y0,
                            std::function<double(double)> f_g,
                            double (&f_s)(double), std::vector<double> &X0,
                            double xi, double xf, double dx) {
    double Nx = fabs(xf - xi) / dx; // Nombre de points de calcul
    std::vector<double> resultat(Nx);

    double x0 = X0[0], x1 = X0[1], x2;
    double y0 = Y0[0], y1 = Y0[1];
    double y2, g0, g1, g2, s0, s1, s2, terme_div;
    resultat[0] = y0;
    resultat[1] = y1;

    int Ind;
    if (x0 < xf) {
        x2 = x1 + dx;
        while (x2 <= xf) {
            g0 = f_g(x0);
            g1 = f_g(x1);
            g2 = f_g(x2);
            s0 = f_s(x0);
            s1 = f_s(x1);
            s2 = f_s(x2);
            terme_div = 1. + dx * dx * g2 / 12.;
            if (terme_div != 0.) {
                y2 = (2 * y1 * (1 - 5 * dx * dx * g1 / 12.) -
                      y0 * (1 + dx * dx * g0 / 12.) +
                      dx * dx * (s2 + 10 * s1 + s0) / 12.) /
                     (terme_div);
            } else {
                std::cerr << "Division par zéro détectée à x2 = " << x2
                          << ". Boucle interrompue." << std::endl;
                break; // Interrompt la boucle en cas de division par zéro
            }

            Ind = (x2 - xi) / dx;
            resultat[Ind] = y2;

            x0 = x1;
            x1 = x2;
            x2 = x1 + dx;
            y0 = y1;
            y1 = y2;
        }

    }

    else if (xf < x0) {
        x2 = x1 - dx;
        while (xf <= x2) {
            g0 = f_g(x0);
            g1 = f_g(x1);
            g2 = f_g(x2);
            s0 = f_s(x0);
            s1 = f_s(x1);
            s2 = f_s(x2);
            terme_div = 1. + dx * dx * g2 / 12.;
            if (terme_div != 0.) {
                y2 = (2 * y1 * (1 - 5 * dx * dx * g1 / 12.) -
                      y0 * (1 + dx * dx * g0 / 12.) +
                      dx * dx * (s2 + 10 * s1 + s0) / 12.) /
                     (terme_div);
            } else {
                std::cerr << "Division par zéro détectée à x2 = " << x2
                          << ". Boucle interrompue." << std::endl;
                break; // Interrompt la boucle en cas de division par zéro
            }

            Ind = -(x2 - xi) / dx;
            resultat[Ind] = y2;

            x0 = x1;
            x1 = x2;
            x2 = x1 - dx;
            y0 = y1;
            y1 = y2;
        }
    }

    return resultat;
}
