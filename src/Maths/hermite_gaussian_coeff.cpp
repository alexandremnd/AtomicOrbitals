#include "Maths/hermite_gaussian_coeff.hpp"
#include "Orbitals/gaussian_primitive.hpp"
#include "Eigen/Dense"

#include <cmath>

bool ConditionIndex(
    int i, int j,
    int t) { // Condition for having the hermite gaussian coefficient non zero
    return (0 <= t && t <= (i + j) && i >= 0 && j >= 0);
}

std::vector<Tensor3D<double>>
HermiteCoefficient(const GaussianPrimitive &orbital1,
                   const GaussianPrimitive &orbital2) {
    double alpha = orbital1.alpha(), beta = orbital2.alpha();
    double p = alpha + beta;
    double mu = alpha * beta / p;
    double AB_proj, PA_proj, PB_proj;
    Eigen::Vector3d A_position = orbital1.position(),
                    B_position = orbital2.position(),
                    AB_position = A_position - B_position;
    Eigen::Vector3d P_position = (alpha * A_position + beta * B_position) / p,
                    PA_position = P_position - A_position,
                    PB_position = P_position - B_position;

    const int i1 = orbital1.x_exponent(), i2 = orbital2.x_exponent();
    const int j1 = orbital1.y_exponent(), j2 = orbital2.y_exponent();
    const int k1 = orbital1.z_exponent(), k2 = orbital2.z_exponent();

    const int iA_loopLimit[] = {i1 + 1, j1 + 1, k1 + 1};
    const int iB_loopLimit[] = {i2 + 1, j2 + 1, k2 + 1};
    const int t_loopLimit[] = {i1 + i2 + 1, j1 + j2 + 1, k1 + k2 + 1};

    std::vector<Tensor3D<double>> E; // Define the list of 3D tensor
    int i, j, t;
    double Eiminus1_j_tminus1, Eiminus1_j_t, Eiminus1_j_tplus1;
    double Ei_jminus1_tminus1, Ei_jminus1_t, Ei_jminus1_tplus1;

    for (int axis = 0; axis <= 2; axis++) {
        Tensor3D<double> E_axis(
            iA_loopLimit[axis], iB_loopLimit[axis],
            t_loopLimit[axis]); // Tensor use to stock hermite gaussian
                                // coefficient on the chosen axis

        AB_proj = AB_position[axis];
        PA_proj = PA_position[axis];
        PB_proj = PB_position[axis];

        E_axis(0, 0, 0) = std::exp(-mu * AB_proj * AB_proj);
        i = 0;
        for (j = 0; j < iB_loopLimit[axis];
             j++) { // Loop for the index of the right gaussian
            for (t = 0; t < t_loopLimit[axis]; t++) {

                if (i != 0 || j != 0 || t != 0) {
                    Ei_jminus1_tminus1 = 0;
                    Ei_jminus1_t = 0;
                    Ei_jminus1_tplus1 = 0;

                    if (ConditionIndex(i, j - 1, t - 1)) {
                        Ei_jminus1_tminus1 = E_axis(i, j - 1, t - 1);
                    }

                    if (ConditionIndex(i, j - 1, t)) {
                        Ei_jminus1_t = E_axis(i, j - 1, t);
                    }

                    if (ConditionIndex(i, j - 1, t + 1)) {
                        Ei_jminus1_tplus1 = E_axis(i, j - 1, t + 1);
                    }

                    E_axis(i, j, t) = Ei_jminus1_tminus1 / (2. * p) +
                                      PB_proj * Ei_jminus1_t +
                                      (t + 1) * Ei_jminus1_tplus1;
                }
            }
        }

        for (i = 1; i < iA_loopLimit[axis];
             i++) { // Loop for the index of the left gaussian
            for (j = 0; j < iB_loopLimit[axis];
                 j++) { // Loop for the index of the right gaussian
                for (t = 0; t < t_loopLimit[axis]; t++) {

                    if (i != 0 || j != 0 || t != 0) {
                        Eiminus1_j_tminus1 = 0;
                        Eiminus1_j_t = 0;
                        Eiminus1_j_tplus1 = 0;

                        if (ConditionIndex(i - 1, j, t - 1)) {
                            Eiminus1_j_tminus1 = E_axis(i - 1, j, t - 1);
                        }

                        if (ConditionIndex(i - 1, j, t)) {
                            Eiminus1_j_t = E_axis(i - 1, j, t);
                        }

                        if (ConditionIndex(i - 1, j, t + 1)) {
                            Eiminus1_j_tplus1 = E_axis(i - 1, j, t + 1);
                        }

                        E_axis(i, j, t) = Eiminus1_j_tminus1 / (2. * p) +
                                          PA_proj * Eiminus1_j_t +
                                          (t + 1) * Eiminus1_j_tplus1;
                    }
                }
            }
        }

        E.push_back(E_axis);
    }
    E.resize(3);
    return E;
}