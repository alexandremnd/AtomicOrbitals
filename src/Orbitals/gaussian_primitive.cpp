#include <boost/math/special_functions/gamma.hpp>
#include <cmath>

#include "Eigen/Dense"
#include "Orbitals/gaussian_primitive.hpp"
#include "Maths/hermite_gaussian_coeff.hpp"
#include "Maths/hermite_integral.hpp"
#include "Maths/tensor3D.hpp"
#include "Maths/tensor4D.hpp"

void GaussianPrimitive::update_normalisation() {
    double fact_x_exponent, fact_y_exponent, fact_z_exponent;
    double fact_2x_exponent, fact_2y_exponent, fact_2z_exponent;

    fact_x_exponent = boost::math::tgamma(m_x_exponent + 1);
    fact_2x_exponent = boost::math::tgamma(2 * m_x_exponent + 1);

    fact_y_exponent = boost::math::tgamma(m_y_exponent + 1);
    fact_2y_exponent = boost::math::tgamma(2 * m_y_exponent + 1);

    fact_z_exponent = boost::math::tgamma(m_z_exponent + 1);
    fact_2z_exponent = boost::math::tgamma(2 * m_z_exponent + 1);

    m_constant =
        std::pow(2 * m_alpha / M_PI, 3. / 4) *
        std::sqrt(
            std::pow(8 * m_alpha, m_x_exponent + m_y_exponent + m_z_exponent) *
            fact_x_exponent * fact_y_exponent * fact_z_exponent /
            fact_2x_exponent / fact_2y_exponent / fact_2z_exponent);
}

void GaussianPrimitive::set_alpha(double alpha) {
    if (alpha <= 0) {
        throw std::invalid_argument(
            "GaussianPrimitive: Alpha must be greater than 0");
    }
    m_alpha = alpha;
}

void GaussianPrimitive::set_x_exponent(int x_exponent) {
    if (x_exponent < 0) {
        throw std::invalid_argument(
            "GaussianPrimitive: x_exponent must be greater or equal than 0");
    }
    m_x_exponent = x_exponent;
}

void GaussianPrimitive::set_y_exponent(int y_exponent) {
    if (y_exponent < 0) {
        throw std::invalid_argument(
            "GaussianPrimitive: z_exponent must be greater or equal than 0");
    }
    m_y_exponent = y_exponent;
}

void GaussianPrimitive::set_z_exponent(int z_exponent) {
    if (z_exponent < 0) {
        throw std::invalid_argument(
            "GaussianPrimitive: z_exponent must be greater or equal than 0");
    }
    m_z_exponent = z_exponent;
}

double GaussianPrimitive::evaluate(double x, double y, double z) const {
    return std::pow(x - m_position(0), m_x_exponent) *
           std::pow(y - m_position(1), m_y_exponent) *
           std::pow(z - m_position(2), m_z_exponent) *
           std::exp(-m_alpha * ((x - m_position(0)) * (x - m_position(0)) +
                                (y - m_position(1)) * (y - m_position(1)) +
                                (z - m_position(2)) * (z - m_position(2))));
}

GaussianPrimitive operator*(const GaussianPrimitive &orbital1,
                            const GaussianPrimitive &orbital2) {

    // Initialize variables with orbitals parameters
    double alpha = orbital1.alpha(), beta = orbital2.alpha();
    double norm1 = orbital1.constant(), norm2 = orbital2.constant();

    Eigen::Vector3d A_position = orbital1.position(),
                    B_position = orbital2.position(),
                    AB_position = A_position - B_position;

    int i = orbital1.x_exponent(), l = orbital2.x_exponent();
    int j = orbital1.y_exponent(), m = orbital2.y_exponent();
    int k = orbital1.z_exponent(), n = orbital2.z_exponent();

    Eigen::Vector3d new_position =
        (alpha * A_position + beta * B_position) / (alpha + beta);

    GaussianPrimitive orbital3 =
        GaussianPrimitive(i + l, j + m, k + n, alpha + beta, new_position);

    orbital3.m_constant = norm1 * norm2 *
                          std::exp(-((alpha * beta) / (alpha + beta)) *
                                   AB_position.norm() * AB_position.norm());

    return orbital3;
}

GaussianPrimitive &GaussianPrimitive::operator=(const GaussianPrimitive &rhs) {
    if (this == &rhs) {
        return *this;
    }

    m_x_exponent = rhs.m_x_exponent;
    m_y_exponent = rhs.m_y_exponent;
    m_z_exponent = rhs.m_z_exponent;
    m_alpha = rhs.m_alpha;
    m_position = rhs.m_position;
    m_constant = rhs.m_constant;

    return *this;
}

// ==============================================================================================
// ========================== Matrix element for Gaussian Primitive
// ==============================================================================================

double overlap_integral(const GaussianPrimitive &orbital1,
                        const GaussianPrimitive &orbital2) {
    std::vector<Tensor3D<double>> Hermite_coeff =
        HermiteCoefficient(orbital1, orbital2);

    double alpha = orbital1.alpha(), beta = orbital2.alpha();

    int i = orbital1.x_exponent(), l = orbital2.x_exponent();
    int j = orbital1.y_exponent(), m = orbital2.y_exponent();
    int k = orbital1.z_exponent(), n = orbital2.z_exponent();

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
    int k = orbital1.z_exponent(), n = orbital2.z_exponent();

    GaussianPrimitive orbital2_plus2(orbital2);
    orbital2_plus2.set_x_exponent(orbital2.x_exponent() + 2);
    orbital2_plus2.set_y_exponent(orbital2.y_exponent() + 2);
    orbital2_plus2.set_z_exponent(orbital2.z_exponent() + 2);

    std::vector<Tensor3D<double>> Hermite_coeff_plus2 =
        HermiteCoefficient(orbital1, orbital2_plus2);
    Tensor3D<double> Herm_x = Hermite_coeff_plus2[0],
                     Herm_y = Hermite_coeff_plus2[1],
                     Herm_z = Hermite_coeff_plus2[2];

    double Til;
    if (l >= 2) {
        Til = 4 * beta * beta * Herm_x(i, l + 2, 0) -
              2 * beta * (2 * l + 1) * Herm_x(i, l, 0) +
              l * (l - 1) * Herm_x(i, l - 2, 0);
    } else {
        Til = 4 * beta * beta * Herm_x(i, l + 2, 0) -
              2 * beta * (2 * l + 1) * Herm_x(i, l, 0);
    }

    double Tjm;
    if (m >= 2) {
        Tjm = 4 * beta * beta * Herm_y(j, m + 2, 0) -
              2 * beta * (2 * m + 1) * Herm_y(j, m, 0) +
              m * (m - 1) * Herm_y(j, m - 2, 0);
    } else {
        Tjm = 4 * beta * beta * Herm_y(j, m + 2, 0) -
              2 * beta * (2 * m + 1) * Herm_y(j, m, 0);
    }

    double Tkn;
    if (n >= 2) {
        Tkn = 4 * beta * beta * Herm_z(k, n + 2, 0) -
              2 * beta * (2 * n + 1) * Herm_z(k, n, 0) +
              n * (n - 1) * Herm_z(k, n - 2, 0);
    } else {
        Tkn = 4 * beta * beta * Herm_z(k, n + 2, 0) -
              2 * beta * (2 * n + 1) * Herm_z(k, n, 0);
    }

    double Eil = Herm_x(i, l, 0), Ejm = Herm_y(j, m, 0), Ekn = Herm_z(k, n, 0);
    double laplacian_value =
        Til * Ejm * Ekn + Eil * Tjm * Ekn + Eil * Ejm * Tkn;

    return std::pow(M_PI / (alpha + beta), 3. / 2) * laplacian_value;
}

double electron_nucleus_integral(const GaussianPrimitive &orbital1,
                                 const GaussianPrimitive &orbital2,
                                 const Eigen::Vector3d &nucleus_position) {
    // Initialize variables with orbitals parameters
    double alpha = orbital1.alpha(), beta = orbital2.alpha();
    double p = alpha + beta;

    Eigen::Vector3d A_position = orbital1.position(),
                    B_position = orbital2.position();

    Eigen::Vector3d P_position = (alpha * A_position + beta * B_position) / p;
    Eigen::Vector3d PNucleus_position = P_position - nucleus_position;

    int i = orbital1.x_exponent(), l = orbital2.x_exponent();
    int j = orbital1.y_exponent(), m = orbital2.y_exponent();
    int k = orbital1.z_exponent(), n = orbital2.z_exponent();

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

    double zeta = p * q / (p + q);

    Eigen::Vector3d position_12_34 =
        orbital12.position() - orbital34.position();
    Tensor4D R_12_34 =
        HermiteIntegral(orbital12, orbital34, zeta, position_12_34);

    int i = orbital12.x_exponent(), l = orbital34.x_exponent();
    int j = orbital12.y_exponent(), m = orbital34.y_exponent();
    int k = orbital12.z_exponent(), n = orbital34.z_exponent();

    int list_x_indices[] = {orbital1.x_exponent(), orbital2.x_exponent(),
                            orbital3.x_exponent(), orbital4.x_exponent()};
    int list_y_indices[] = {orbital1.y_exponent(), orbital2.y_exponent(),
                            orbital3.y_exponent(), orbital4.y_exponent()};
    int list_z_indices[] = {orbital1.z_exponent(), orbital2.z_exponent(),
                            orbital3.z_exponent(), orbital4.z_exponent()};

    double E12x, E12y, E12z, E34x, E34y, E34z, R0_12_34;
    double electron_electron_value = 0;
    double E_product;
    for (int t = 0; t <= i; t++) {
        E12x = E12[0](list_x_indices[0], list_x_indices[1], t);
        for (int u = 0; u <= j; u++) {
            E12y = E12[1](list_y_indices[0], list_y_indices[1], u);
            for (int v = 0; v <= k; v++) {
                E12z = E12[2](list_z_indices[0], list_z_indices[1], v);
                for (int tau = 0; tau <= l; tau++) {
                    E34x = E34[0](list_x_indices[2], list_x_indices[3], tau);
                    for (int mu = 0; mu <= m; mu++) {
                        E34y = E34[1](list_y_indices[2], list_y_indices[3], mu);
                        for (int nu = 0; nu <= n; nu++) {
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

                            // std::cout << R0_12_34 << std::endl;
                            electron_electron_value +=
                                std::pow(-1, tau + mu + nu) * E_product *
                                R0_12_34;
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