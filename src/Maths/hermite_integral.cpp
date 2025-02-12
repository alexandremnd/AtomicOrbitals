#include "Maths/hermite_integral.hpp"
#include "Eigen/Dense"
#include "Maths/tensor4D.hpp"
#include "Maths/the_boys.hpp"

#include <cmath>
#include <vector>

bool ConditionRecurrence(int n, int t, int u, int v) {
    return (n >= 0 && t >= 0 && u >= 0 && v >= 0);
}

Tensor4D<double> HermiteIntegral(const GaussianPrimitive &orbital1,
                                 const GaussianPrimitive &orbital2, double p,
                                 Eigen::Vector3d position) {

    int i = orbital1.x_exponent(), l = orbital2.x_exponent();
    int j = orbital1.y_exponent(), m = orbital2.y_exponent();
    int k = orbital1.z_exponent(), n = orbital2.z_exponent();
    int ijkMax = i + l + j + m + k + n + 1;
    int loop_limit[] = {i + l + 1, j + m + 1, k + n + 1};

    Tensor4D<double> Rntuv(ijkMax, loop_limit[0], loop_limit[1],
                           loop_limit[2]); // Define the 4 dimension tensor
    double evaluate_term =
        p * position.norm() *
        position.norm(); // zeta is the multiplicating term with the distance,
                         // equal to alpha + beta for electron-nucleus but not
                         // for electron-electron
    std::vector<double> Fn_values = TheBoys_Recurrence(ijkMax, evaluate_term);

    int n_ind, t, u, v;
    for (n_ind = 0; n_ind < ijkMax; n_ind++) {
        Rntuv(n_ind, 0, 0, 0) =
            std::pow(-2 * p, n_ind) *
            Fn_values[n_ind]; // Initialize the values of Rn000 with The Boys
                              // functions
    }

    double Rminus2, Rminus1;
    int maxtuv;
    for (t = 0; t < loop_limit[0]; t++) {
        for (u = 0; u < loop_limit[1]; u++) {
            for (v = 0; v < loop_limit[2]; v++) {
                for (n_ind = 0; n_ind < t + u + v + 1; n_ind++) {

                    maxtuv = std::max(t, std::max(u, v));
                    if (t == 0 && u == 0 && v == 0) {
                        continue;
                    }

                    else if (t == maxtuv) {
                        if (ConditionRecurrence(n_ind + 1, t - 2, u, v)) {
                            Rminus2 = Rntuv(n_ind + 1, t - 2, u, v);
                        } else {
                            Rminus2 = 0;
                        }

                        if (ConditionRecurrence(n_ind + 1, t - 1, u, v)) {
                            Rminus1 = Rntuv(n_ind + 1, t - 1, u, v);
                        } else {
                            Rminus1 = 0;
                        }

                        Rntuv(n_ind, t, u, v) =
                            (t - 1) * Rminus2 + position.x() * Rminus1;
                    }

                    else if (u == maxtuv) {
                        if (ConditionRecurrence(n_ind + 1, t, u - 2, v)) {
                            Rminus2 = Rntuv(n_ind + 1, t, u - 2, v);
                        } else {
                            Rminus2 = 0;
                        }

                        if (ConditionRecurrence(n_ind + 1, t, u - 1, v)) {
                            Rminus1 = Rntuv(n_ind + 1, t, u - 1, v);
                        } else {
                            Rminus1 = 0;
                        }

                        Rntuv(n_ind, t, u, v) =
                            (u - 1) * Rminus2 + position.y() * Rminus1;
                    }

                    else if (v == maxtuv) {
                        if (ConditionRecurrence(n_ind + 1, t, u, v - 2)) {
                            Rminus2 = Rntuv(n_ind + 1, t, u, v - 2);
                        } else {
                            Rminus2 = 0;
                        }

                        if (ConditionRecurrence(n_ind + 1, t, u, v - 1)) {
                            Rminus1 = Rntuv(n_ind + 1, t, u, v - 1);
                        } else {
                            Rminus1 = 0;
                        }

                        Rntuv(n_ind, t, u, v) =
                            (v - 1) * Rminus2 + position.z() * Rminus1;
                    }
                }
            }
        }
    }
    // Rntuv.dim();
    return Rntuv;
}