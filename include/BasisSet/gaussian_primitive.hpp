#pragma once

#include <boost/math/special_functions/gamma.hpp>
#include <cmath>

#include "BasisSet/orbital.hpp"
#include "Eigen/Dense"

class GaussianPrimitive : public Orbital {
  public:
    GaussianPrimitive() {};

    GaussianPrimitive(int x_exponent, int y_exponent, int z_exponent,
                      double alpha)
        : m_x_exponent(x_exponent), m_y_exponent(y_exponent),
          m_z_exponent(z_exponent), m_alpha(alpha) {
        m_position = Eigen::Vector3d::Zero();

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
        m_normalization =
            std::pow(2 * m_alpha / M_PI, 3. / 4) *
            std::sqrt(std::pow(8 * m_alpha,
                               m_x_exponent + m_y_exponent + m_z_exponent) *
                      fact_x_exponent * fact_y_exponent * fact_z_exponent /
                      fact_2x_exponent / fact_2y_exponent / fact_2z_exponent);
    }

    GaussianPrimitive(int x_exponent, int y_exponent, int z_exponent,
                      double alpha, Eigen::Vector3d position)
        : GaussianPrimitive(x_exponent, y_exponent, z_exponent, alpha) {
        m_position = position;
    }

    GaussianPrimitive(GaussianPrimitive const &other)
        : m_x_exponent(other.m_x_exponent), m_y_exponent(other.m_y_exponent),
          m_z_exponent(other.m_z_exponent), m_alpha(other.m_alpha) {
        m_position = other.m_position;
        m_normalization = other.m_normalization;
    }

    GaussianPrimitive(GaussianPrimitive &&other)
        : m_x_exponent(other.m_x_exponent), m_y_exponent(other.m_y_exponent),
          m_z_exponent(other.m_z_exponent), m_alpha(other.m_alpha) {
        m_position = std::move(other.m_position);
        m_normalization = other.m_normalization;
    }

    double evaluate(double x, double y, double z) const;
    double evaluate(double r) const;

    // Setters & Getters
    void set_normalization(double normalisation) {
        m_normalization = normalisation;
    }
    void set_alpha(double alpha) { m_alpha = alpha; }
    void set_x_exponent(int x_exponent) { m_x_exponent = x_exponent; }
    void set_y_exponent(int y_exponent) { m_y_exponent = y_exponent; }
    void set_z_exponent(int z_exponent) { m_z_exponent = z_exponent; }
    void set_position(Eigen::Vector3d position) { m_position = position; }

    int x_exponent() const { return m_x_exponent; }
    int y_exponent() const { return m_y_exponent; }
    int z_exponent() const { return m_z_exponent; }
    double alpha() const { return m_alpha; }
    double normalization() const { return m_normalization; }
    Eigen::Vector3d position() const { return m_position; }

    friend GaussianPrimitive operator*(const GaussianPrimitive &orbital1,
                                       const GaussianPrimitive &orbital2);

  private:
    int m_x_exponent;
    int m_y_exponent;
    int m_z_exponent;
    double m_alpha, m_normalization;
};

double overlap_integral(const GaussianPrimitive &, const GaussianPrimitive &);
double laplacian_integral(const GaussianPrimitive &, const GaussianPrimitive &);
double electron_nucleus_integral(const GaussianPrimitive &,
                                 const GaussianPrimitive &,
                                 const Eigen::Vector3d &);
double electron_electron_integral(const GaussianPrimitive &,
                                  const GaussianPrimitive &,
                                  const GaussianPrimitive &,
                                  const GaussianPrimitive &);