#pragma once

#include <vector>

#include "Orbitals/gaussian_primitive.hpp"

/**
 * @brief Linear combination of primitive basis functions.
 *
    \f[
        \Psi = \sum_i c_i \phi_i
    \f]
    where \f$ c_i \f$ are the coefficients and \f$ \phi_i \f$ are the primitive
 basis functions.
 * @note Coefficients needs to be normalized. Otherwise, you should add
 * primitive and call update_normalization after adding all primitives.
 */
class ContractedGaussian final : public Orbital {
  public:
    ContractedGaussian() = default;

    /**
     * @param size Number of primitives for the linear combination.
     */
    ContractedGaussian(size_t size)
        : m_coefficients(size), m_primitives(size) {}

    /**
     * @param coefficients Weights for each primitive in the linear combination.
     * @param primitives Primitives basis functions.
     * @throws std::length_error if the number of coefficients and primitives
     * are not equal.
     */
    ContractedGaussian(std::vector<double> &coefficients,
                       std::vector<GaussianPrimitive> &primitives);

    ContractedGaussian(const ContractedGaussian &other);

    ContractedGaussian(const ContractedGaussian &&other);

    /**
     * @brief Allocates memory for primitives and coefficients.
     *
     * @param size Number of primitives for the linear combination.
     */
    void reserve(std::size_t size);

    /**
     * @brief Adds a primitive to the linear combination.
     *
     * @param coefficient Weight of the primitive in the linear combination.
     * @param primitive Primitive to add.
     */
    void add_primitive(double coefficient, const GaussianPrimitive &primitive);

    /**
     * @brief Adds a primitive to the linear combination.
     *
     * @tparam Args
     * @param coefficient Weight of the primitive in the linear combination.
     * @param args Arguments to pass to the GaussianPrimitive constructor.
     */
    template <typename... Args>
    void add_primitive(double coefficient, Args &&...args) {
        m_coefficients.push_back(coefficient);
        m_primitives.emplace_back(std::forward<Args>(args)...);
    }

    void update_normalization();

    size_t size() const;

    const std::vector<GaussianPrimitive> &get_primitives() const;
    const GaussianPrimitive &get_primitive(int i) const;
    double get_coefficient(int i) const;

    void set_position(Eigen::Vector3d position) override {
        m_position = position;

        for (auto &primitive : m_primitives) {
            primitive.set_position(position);
        }
    }

    friend std::ostream &operator<<(std::ostream &os,
                                    const ContractedGaussian &orbital) {
        os << "Nucleus position: " << orbital.position().transpose()
           << std::endl;
        for (size_t i = 0; i < orbital.size(); i++) {
            os << "\t * Primitive" << i << ": "
               << orbital.get_coefficient(i) *
                      orbital.get_primitive(i).constant()
               << orbital.get_primitive(i) << std::endl;
        }
        return os;
    }

  private:
    std::vector<double> m_coefficients;
    std::vector<GaussianPrimitive> m_primitives;
};

typedef ContractedGaussian CGTO;

/**
 * @brief Computes the overlap integral between two contracted orbitals <o1|o2>
 *
 * @param o1 First contracted orbital
 * @param o2 Second contracted orbital
 * @return double <o1|o2>
 */

double overlap_integral(const ContractedGaussian &o1,
                        const ContractedGaussian &o2);

/**
 * @brief Computes the laplacian integral between two contracted orbitals
 *
 * @param o1 First contracted orbital
 * @param o2 Second contracted orbital
 * @return double <o1|nabla^2|o2>
 */

double laplacian_integral(const ContractedGaussian &o1,
                          const ContractedGaussian &o2);

/**
 * @brief Computes the 1/|r-pos| integral between two contracted orbitals
 *
 * @param o1 First contracted orbital
 * @param o2 Second contracted orbital
 * @param pos Position of the considered nucleus
 * @return double <o1|1/|r-pos| |o2>
 */

double electron_nucleus_integral(const ContractedGaussian &o1,
                                 const ContractedGaussian &o2,
                                 const Eigen::Vector3d &pos);

/**
 * @brief Computes the electron-electron integral between four contractd
 * orbitals
 *
 * @param o1 First contracted orbital
 * @param o2 Second contracted orbital sharing coordinates of o1
 * @param o3 Third contracted orbital
 * @param o4 Fourth contracted orbital sharing coordinates of o3
 * @return double <o1 o3 | 1/|r1-r2| | o2 o4>
 */

double electron_electron_integral(const ContractedGaussian &o1,
                                  const ContractedGaussian &o2,
                                  const ContractedGaussian &o3,
                                  const ContractedGaussian &o4);
