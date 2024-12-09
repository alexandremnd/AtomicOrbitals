#pragma once

#include "Eigen/Dense"

class GaussianPrimitive {
public:
    GaussianPrimitive(int x_exponent, int y_exponent, int z_exponent, double alpha) : m_x_exponent(x_exponent), m_y_exponent(y_exponent), m_z_exponent(z_exponent), m_alpha(alpha) {
        m_normalization = 1.0;
    }

    GaussianPrimitive(int x_exponent, int y_exponent, int z_exponent, double alpha, Eigen::Vector3d position) : GaussianPrimitive(x_exponent, y_exponent, z_exponent, alpha) {
        m_position = position;
    }

    GaussianPrimitive(GaussianPrimitive const &other) : m_x_exponent(other.m_x_exponent), m_y_exponent(other.m_y_exponent), m_z_exponent(other.m_z_exponent), m_alpha(other.m_alpha), m_normalization(other.m_normalization), m_position(other.m_position) {}

    GaussianPrimitive(GaussianPrimitive &&other) : m_x_exponent(other.m_x_exponent), m_y_exponent(other.m_y_exponent), m_z_exponent(other.m_z_exponent), m_alpha(other.m_alpha), m_normalization(other.m_normalization) {
        m_position = std::move(other.m_position);
    }

    double evaluate(double x, double y, double z) const;
    double evaluate(double r) const;

    // Getters
    void set_x_exponent(int x_exponent) { m_x_exponent = x_exponent; }
    void set_y_exponent(int y_exponent) { m_y_exponent = y_exponent; }
    void set_z_exponent(int z_exponent) { m_z_exponent = z_exponent; }

    int x_exponent() const { return m_x_exponent; }
    int y_exponent() const { return m_y_exponent; }
    int z_exponent() const { return m_z_exponent; }
    double alpha() const { return m_alpha; }
    double normalization() const { return m_normalization; }

private:
    int m_x_exponent;
    int m_y_exponent;
    int m_z_exponent;
    double m_alpha, m_normalization;
    Eigen::Vector3d m_position;
};