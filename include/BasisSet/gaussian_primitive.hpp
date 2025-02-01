#pragma once

#include "Eigen/Dense"
#include "BasisSet/orbital.hpp"

class GaussianPrimitive final : public Orbital {
public:
    GaussianPrimitive(int x_exponent, int y_exponent, int z_exponent, double alpha, Eigen::Vector3d position) :
    m_x_exponent(x_exponent), m_y_exponent(y_exponent), m_z_exponent(z_exponent), m_alpha(alpha) {
        m_position = position;
        m_normalization_constant = 1.0;
    }

    GaussianPrimitive() : GaussianPrimitive(0, 0, 0, 0.0, Eigen::Vector3d{0, 0, 0}) {}

    GaussianPrimitive(int x_exponent, int y_exponent, int z_exponent, double alpha) : GaussianPrimitive(x_exponent, y_exponent, z_exponent, alpha, Eigen::Vector3d{0, 0, 0}) {}

    GaussianPrimitive(GaussianPrimitive const &other) = default;

    double evaluate(double x, double y, double z) const;

    void set_x_exponent(int x_exponent);
    void set_y_exponent(int y_exponent);
    void set_z_exponent(int z_exponent);

    // Getters
    inline int x_exponent() const { return m_x_exponent; }
    inline int y_exponent() const { return m_y_exponent; }
    inline int z_exponent() const { return m_z_exponent; }
    inline double alpha() const { return m_alpha; }

private:
    int m_x_exponent;
    int m_y_exponent;
    int m_z_exponent;
    double m_alpha;
};