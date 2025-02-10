#include "BasisSet/gaussian_primitive.hpp"
#include "Maths/hermite_gaussian_coeff.hpp"
#include "Maths/hermite_integral.hpp"
#include "Maths/tensor3D.hpp"
#include "Maths/tensor4D.hpp"

double GaussianPrimitive::evaluate(double x, double y, double z)
    const { // Return x^i * y^j * z^k * exp(-alpha*[(x-xA)² + (y-yA)² +
            // (z-zA)²])
    return std::pow(x - m_position(0), m_x_exponent) *
           std::pow(y - m_position(1), m_y_exponent) *
           std::pow(z - m_position(2), m_z_exponent) *
           std::exp(-m_alpha * ((x - m_position(0)) * (x - m_position(0)) +
                                (y - m_position(1)) * (y - m_position(1)) +
                                (z - m_position(2)) * (z - m_position(2))));
}

GaussianPrimitive operator*(
    const GaussianPrimitive &orbital1,
    const GaussianPrimitive &orbital2) { // return the product of two gaussian
                                         // primitives with the same position
    // Initialize variables with orbitals parameters
    double alpha = orbital1.alpha(), beta = orbital2.alpha();
    double norm1 = orbital1.normalization(), norm2 = orbital2.normalization();

    Eigen::Vector3d A_position = orbital1.position(),
                    B_position = orbital2.position();

    int i = orbital1.x_exponent(), l = orbital2.x_exponent();
    int j = orbital1.y_exponent(), m = orbital2.y_exponent();
    int k = orbital1.y_exponent(), n = orbital2.y_exponent();

    GaussianPrimitive orbital3;
    if (A_position == B_position) {
        orbital3.set_position(A_position);
    }
    orbital3.set_normalization(norm1 * norm2);
    orbital3.set_alpha(alpha + beta);
    orbital3.set_x_exponent(i + l);
    orbital3.set_y_exponent(j + m);
    orbital3.set_z_exponent(k + n);

    return orbital3;
}

// ===============================================================================================
// ============================ Matrix element for Gaussian Primitive
// ============================
// ===============================================================================================

double overlap_integral(const GaussianPrimitive &orbital1,
                        const GaussianPrimitive &orbital2) {
    std::vector<Tensor3D<double>> Hermite_coeff =
        HermiteCoefficient(orbital1, orbital2);

    double alpha = orbital1.alpha(), beta = orbital2.alpha();

    int i = orbital1.x_exponent(), l = orbital2.x_exponent();
    int j = orbital1.y_exponent(), m = orbital2.y_exponent();
    int k = orbital1.y_exponent(), n = orbital2.y_exponent();

    double Eil0 = Hermite_coeff[0](i, l, 0), Ejm0 = Hermite_coeff[1](j, m, 0),
           Ekn0 = Hermite_coeff[2](k, n, 0);
    double overlap_value =
        std::pow(M_PI / (alpha + beta), 3. / 2) * Eil0 * Ejm0 * Ekn0;
    return overlap_value;
}

double laplacian_integral(const GaussianPrimitive &orbital1,
                          const GaussianPrimitive &orbital2) {
    double alpha = orbital1.alpha(), beta = orbital2.alpha();

    int i = orbital1.x_exponent(), l = orbital2.x_exponent();
    int j = orbital1.y_exponent(), m = orbital2.y_exponent();
    int k = orbital1.y_exponent(), n = orbital2.y_exponent();

    std::vector<Tensor3D<double>> Hermite_coeff =
        HermiteCoefficient(orbital1, orbital2);
    Tensor3D Coeff_xaxis = Hermite_coeff[0], Coeff_yaxis = Hermite_coeff[1],
             Coeff_zaxis = Hermite_coeff[2];

    double Eil0 = Coeff_xaxis(i, l, 0), Ejm0 = Coeff_xaxis(j, m, 0),
           Ekn0 = Coeff_xaxis(k, n, 0),
           Til = 4 * beta * beta * Coeff_xaxis(i, l + 2, 0) -
                 2 * beta * (2 * l + 1) * Eil0 +
                 l * (l - 1) * Coeff_xaxis(i, l - 2, 0),
           Tjm = 4 * beta * beta * Coeff_yaxis(j, m + 2, 0) -
                 2 * beta * (2 * m + 1) * Ejm0 +
                 m * (m - 1) * Coeff_yaxis(j, m - 2, 0),
           Tkn = 4 * beta * beta * Coeff_zaxis(k, n + 2, 0) -
                 2 * beta * (2 * k + 1) * Ekn0 +
                 k * (k - 1) * Coeff_zaxis(k, n - 2, 0),

           laplacian_value =
               -(1. / 2) *
               (Til * Ejm0 * Ekn0 + Eil0 * Tjm * Ekn0 + Eil0 * Ejm0 * Tkn) *
               std::pow(M_PI / (alpha + beta), 3. / 2);

    return laplacian_value;
}

double electron_nucleus_integral(const GaussianPrimitive &orbital1,
                                 const GaussianPrimitive &orbital2,
                                 const Eigen::Vector3d &nucleus_position) {
    // Initialize variables with orbitals parameters
    double alpha = orbital1.alpha(), beta = orbital2.alpha();
    double p = alpha + beta;

    Eigen::Vector3d A_position = orbital1.position(),
                    B_position = orbital2.position(),
                    AB_position = A_position - B_position;
    Eigen::Vector3d P_position = (alpha * A_position + beta * B_position) / p;
    Eigen::Vector3d PNucleus_position = P_position - nucleus_position;

    int i = orbital1.x_exponent(), l = orbital2.x_exponent();
    int j = orbital1.y_exponent(), m = orbital2.y_exponent();
    int k = orbital1.y_exponent(), n = orbital2.y_exponent();
    int ijkMax = i + l + j + m + k + n;

    // Initialize variables for Hermite coefficients and Hermite integrals
    std::vector<Tensor3D<double>> Hermite_coeff =
        HermiteCoefficient(orbital1, orbital2);
    Tensor3D E_xaxis = Hermite_coeff[0], E_yaxis = Hermite_coeff[1],
             E_zaxis = Hermite_coeff[2];
    Tensor4D Rntuv =
        HermiteIntegral(orbital1, orbital2, p,
                        PNucleus_position); // Define the 4 dimension tensor

    double Eilt, Ejmu, Eknv, R0tuv;
    double electron_nucleus_value = 0;
    for (int t = 0; t <= i + l; t++) {
        for (int u = 0; u <= j + m; u++) {
            for (int v = 0; v <= k + n; v++) {
                Eilt = E_xaxis(i, l, t);
                Ejmu = E_yaxis(j, m, u);
                Eknv = E_zaxis(k, n, v);
                R0tuv = Rntuv(0, t, u, v);
                electron_nucleus_value += Eilt * Ejmu * Eknv * R0tuv;
            }
        }
    }

    electron_nucleus_value *= 2 * M_PI / p;
    return electron_nucleus_value;
}

double electron_electron_integral(const GaussianPrimitive &orbital1,
                                  const GaussianPrimitive &orbital2,
                                  const GaussianPrimitive &orbital3,
                                  const GaussianPrimitive &orbital4) {
    std::vector<Tensor3D<double>> E12 = HermiteCoefficient(orbital1, orbital2),
                                  E34 = HermiteCoefficient(orbital3, orbital4);
    GaussianPrimitive orbital12 = orbital1 * orbital2,
                      orbital34 = orbital3 * orbital4;
    double p = orbital12.alpha(), q = orbital34.alpha();

    int zeta = p * q / (p + q);
    Eigen::Vector3d position_12_34 =
        orbital12.position() - orbital34.position();
    Tensor4D R_12_34 =
        HermiteIntegral(orbital12, orbital34, zeta, position_12_34);

    int i = orbital12.x_exponent(), l = orbital34.x_exponent();
    int j = orbital12.y_exponent(), m = orbital34.y_exponent();
    int k = orbital12.y_exponent(), n = orbital34.y_exponent();

    int list_x_indices[] = {orbital1.x_exponent(), orbital2.x_exponent(),
                            orbital3.x_exponent(), orbital4.x_exponent()};
    int list_y_indices[] = {orbital1.y_exponent(), orbital2.y_exponent(),
                            orbital3.y_exponent(), orbital4.y_exponent()};
    int list_z_indices[] = {orbital1.z_exponent(), orbital2.z_exponent(),
                            orbital3.z_exponent(), orbital4.z_exponent()};

    double E12x, E12y, E12z, E34x, E34y, E34z, R0_12_34;
    double electron_electron_value = 0;
    double E_product;

    for (int t; t <= i; t++) {
        E12x = E12[0](list_x_indices[0], list_x_indices[1], t);
        for (int u; u <= j; u++) {
            E12y = E12[1](list_y_indices[0], list_y_indices[1], u);
            for (int v; v <= k; v++) {
                E12z = E12[2](list_z_indices[0], list_z_indices[1], v);
                for (int tau; tau <= l; tau++) {
                    E34x = E34[0](list_x_indices[2], list_x_indices[3], tau);
                    for (int mu; mu <= m; mu++) {
                        E34y = E34[1](list_y_indices[2], list_y_indices[3], mu);
                        for (int nu; nu <= n; nu++) {
                            E34z = E34[2](list_z_indices[2], list_z_indices[3],
                                          nu);

                            E_product = 1;
                            E_product *= E12x;
                            E_product *= E12y;
                            E_product *= E12z;

                            E_product *= E34x;
                            E_product *= E34y;
                            E_product *= E34z;

                            R0_12_34 = R_12_34(0, t + tau, u + mu, v + nu);

                            electron_electron_value += E_product * R0_12_34;
                        }
                    }
                }
            }
        }
    }

    electron_electron_value *=
        2 * std::pow(M_PI, 5. / 2) / (p * q * std::sqrt(p + q));
    return electron_electron_value;
}