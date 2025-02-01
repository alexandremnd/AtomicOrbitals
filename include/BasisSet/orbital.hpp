#pragma once

#include "Eigen/Dense"

class Orbital {
    public:
        ~Orbital() = default;

        double normalization() const { return m_normalization_constant; };
        Eigen::Vector3d position() const { return m_position; };

        virtual void set_position(Eigen::Vector3d position) { m_position = position; };

    protected:
        double m_normalization_constant = 0;
        Eigen::Vector3d m_position = Eigen::Vector3d{0, 0, 0};
};