#pragma once

#include <boost/math/special_functions/gamma.hpp>
#include <cmath>

#include "Orbitals/orbital.hpp"
#include "Eigen/Dense"

/**
 * @brief Gaussian primitive orbital.
 * Orbital defined as:
 * \f[
 *      \phi(x, y, z) = N \cdot (x-x_A)^i (y-y_A)^j (z-z_A)^k \cdot e^{-\alpha
 *      ||\mathbf{r} - \mathbf{R_A}||^2}
 * \f]
 * where \f$N\f$ is the normalization constant, \f$\alpha\f$ is the decay rate,
 * \f$i, j, k\f$ are the exponents in the x, y, z directions respectively, and
 * \f$\mathbf{R_A}\f$ is the position of the nucleus.
 *
 */
class GaussianPrimitive : public Orbital {
  public:
    GaussianPrimitive(){};

    GaussianPrimitive(int x_exponent, int y_exponent, int z_exponent,
                      double alpha, Eigen::Vector3d position)
        : m_x_exponent(x_exponent), m_y_exponent(y_exponent),
          m_z_exponent(z_exponent), m_alpha(alpha) {
        m_position = position;
        update_normalisation();
    }

    GaussianPrimitive(int x_exponent, int y_exponent, int z_exponent,
                      double alpha)
        : GaussianPrimitive(x_exponent, y_exponent, z_exponent, alpha,
                            Eigen::Vector3d{0., 0., 0.}) {}

    GaussianPrimitive(GaussianPrimitive const &other)
        : m_x_exponent(other.m_x_exponent), m_y_exponent(other.m_y_exponent),
          m_z_exponent(other.m_z_exponent), m_alpha(other.m_alpha) {
        m_position = other.m_position;
        m_normalization_constant = other.m_normalization_constant;
    }

    GaussianPrimitive(GaussianPrimitive &&other)
        : m_x_exponent(other.m_x_exponent), m_y_exponent(other.m_y_exponent),
          m_z_exponent(other.m_z_exponent), m_alpha(other.m_alpha) {
        m_position = other.m_position;
        m_normalization_constant = other.m_normalization_constant;
    }

    /**
     * @brief Evaluates the orbital at a given position.
     *
     * @param x x-coordinate
     * @param y y-coordinate
     * @param z z-coordinate
     * @return double Phi(x, y, z)
     */
    double evaluate(double x, double y, double z) const;

    void update_normalisation() {
        double fact_x_exponent, fact_y_exponent, fact_z_exponent;
        double fact_2x_exponent, fact_2y_exponent, fact_2z_exponent;

        if (m_x_exponent == 0) {
            fact_x_exponent = 1.;
            fact_2x_exponent = 1.;
        } else {
            fact_x_exponent = boost::math::tgamma(m_x_exponent + 1);
            fact_2x_exponent = boost::math::tgamma(2 * m_x_exponent + 1);
        }

        if (m_y_exponent == 0) {
            fact_y_exponent = 1.;
            fact_2y_exponent = 1.;
        } else {
            fact_y_exponent = boost::math::tgamma(m_y_exponent + 1);
            fact_2y_exponent = boost::math::tgamma(2 * m_y_exponent + 1);
        }

        if (m_z_exponent == 0) {
            fact_z_exponent = 1.;
            fact_2z_exponent = 1.;
        } else {
            fact_z_exponent = boost::math::tgamma(m_z_exponent + 1);
            fact_2z_exponent = boost::math::tgamma(2 * m_z_exponent + 1);
        }
        m_normalization_constant =
            std::pow(2 * m_alpha / M_PI, 3. / 4) *
            std::sqrt(std::pow(8 * m_alpha,
                               m_x_exponent + m_y_exponent + m_z_exponent) *
                      fact_x_exponent * fact_y_exponent * fact_z_exponent /
                      fact_2x_exponent / fact_2y_exponent / fact_2z_exponent);
    }

    // Setters & Getters
    /**
     * @throw std::invalid_argument if alpha <= 0
     */
    void set_alpha(double alpha);

    /**
     * @throw std::invalid_argument if x_exponent < 0
     */
    void set_x_exponent(int x_exponent);

    /**
     * @throw std::invalid_argument if y_exponent < 0
     */
    void set_y_exponent(int y_exponent);

    /**
     * @throw std::invalid_argument if z_exponent < 0
     */
    void set_z_exponent(int z_exponent);

    int x_exponent() const { return m_x_exponent; }
    int y_exponent() const { return m_y_exponent; }
    int z_exponent() const { return m_z_exponent; }
    double alpha() const { return m_alpha; }

    friend GaussianPrimitive operator*(const GaussianPrimitive &orbital1,
                                       const GaussianPrimitive &orbital2);

    GaussianPrimitive &operator=(const GaussianPrimitive &rhs);

  private:
    int m_x_exponent;
    int m_y_exponent;
    int m_z_exponent;
    double m_alpha;
};

typedef GaussianPrimitive GTO;

double overlap_integral(const GaussianPrimitive &, const GaussianPrimitive &);
double laplacian_integral(const GaussianPrimitive &, const GaussianPrimitive &);
double electron_nucleus_integral(const GaussianPrimitive &,
                                 const GaussianPrimitive &,
                                 const Eigen::Vector3d &);
double electron_electron_integral(const GaussianPrimitive &,
                                  const GaussianPrimitive &,
                                  const GaussianPrimitive &,
                                  const GaussianPrimitive &);