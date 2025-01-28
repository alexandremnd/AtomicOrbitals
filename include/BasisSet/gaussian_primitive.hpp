#pragma once

#include "Eigen/Dense"
#include "BasisSet/orbital.hpp"

class GaussianPrimitive final : public Orbital {
public:
    GaussianPrimitive() : m_x_exponent(0), m_y_exponent(0), m_z_exponent(0), m_alpha(0.0), m_normalization(0.0) {}

    GaussianPrimitive(int x_exponent, int y_exponent, int z_exponent, double alpha) : m_x_exponent(x_exponent), m_y_exponent(y_exponent), m_z_exponent(z_exponent), m_alpha(alpha) {
        m_normalization_constant = 1.0;
    }

    GaussianPrimitive(int x_exponent, int y_exponent, int z_exponent, double alpha, Eigen::Vector3d position) : GaussianPrimitive(x_exponent, y_exponent, z_exponent, alpha) {
        m_position = position;
    }

    GaussianPrimitive(GaussianPrimitive const &other) : m_x_exponent(other.m_x_exponent), m_y_exponent(other.m_y_exponent), m_z_exponent(other.m_z_exponent), m_alpha(other.m_alpha), m_normalization(other.m_normalization), m_position(other.m_position) {}

    GaussianPrimitive(GaussianPrimitive &&other) : m_x_exponent(other.m_x_exponent), m_y_exponent(other.m_y_exponent), m_z_exponent(other.m_z_exponent), m_alpha(other.m_alpha), m_normalization(other.m_normalization) {
        m_position = std::move(other.m_position);
    }

    double evaluate(double x, double y, double z) const;

    // Getters
    void set_x_exponent(int x_exponent);
    void set_y_exponent(int y_exponent);
    void set_z_exponent(int z_exponent);

    inline int x_exponent() const { return m_x_exponent; }
    inline int y_exponent() const { return m_y_exponent; }
    inline int z_exponent() const { return m_z_exponent; }
    inline double alpha() const { return m_alpha; }

private:
    int m_x_exponent;
    int m_y_exponent;
    int m_z_exponent;
    double m_alpha;
    Eigen::Vector3d m_position;
};