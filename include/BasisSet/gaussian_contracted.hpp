#pragma once

#include <vector>
#include "BasisSet/gaussian_primitive.hpp"
#include "BasisSet/orbital.hpp"

class ContractedGaussian  final : public Orbital {
public:
    ContractedGaussian(const std::vector<double>& coefficients,
                       const std::vector<GaussianPrimitive>& primitives)
        : m_coefficients(coefficients), m_primitives(primitives) {}

    ContractedGaussian(const int size) : m_coefficients(size), m_primitives(size) {}

    ContractedGaussian() = default;

    inline double evaluate(double x, double y, double z) const;

    /**
     * @brief Add a primitive to the contracted gaussian orbital.
     *
     * @param coefficient Coefficient of the primitive in the linear combination.
     * @param decay Exponential decay rate of the primitive.
     * @param x_exponent x exponent of the primitive.
     * @param y_exponent y exponent of the primitive.
     * @param z_exponent z exponent of the primitive.
     */
    void add_primitive(double coefficient, double decay, int x_exponent, int y_exponent, int z_exponent, Eigen::Vector3d position = {0, 0, 0}) {
        m_coefficients.push_back(coefficient);
        m_primitives.push_back(GaussianPrimitive(x_exponent, y_exponent, z_exponent, decay, position));
    }

    /**
     * @brief Get the size of the contracted gaussian (number of primitives).
     *
     * @return int Number of primitives in the contracted gaussian.
     */
    inline int get_gaussian_count() const {
        return m_primitives.size();
    }

    /**
     * @brief Get the coefficient of a primitive in the contracted gaussian.
     *
     * @param i i-th primitive.
     * @return double Coefficient of the i-th primitive.
     */
    inline double get_coefficient(int i) const {
        return m_coefficients[i];
    }

    /**
     * @brief Get the primitive at index i.
     *
     * @param i i-th primitive.
     * @return const GaussianPrimitive& Reference to the i-th primitive.
     */
    inline const GaussianPrimitive& get_primitive(int i) const {
        return m_primitives[i];
    }

private:
    std::vector<double> m_coefficients;
    std::vector<GaussianPrimitive> m_primitives;
};