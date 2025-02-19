#pragma once

#include <vector>

#include "Orbitals/gaussian_primitive.hpp"
#include "Orbitals/slater_primitive.hpp"
#include "concepts.hpp"

/**
 * @brief Linear combination of primitive basis functions.
 *
    \f[
        \Psi = \sum_i c_i \phi_i
    \f]
    where \f$ c_i \f$ are the coefficients and \f$ \phi_i \f$ are the primitive
 basis functions.
 * @note Coefficients do not need to be normalized, the normalization constant
 is calculated automatically when adding a primitive.
 *
 * @tparam PrimitiveType Type of primitive basis functions (e.g
 GaussianPrimitive, SlaterPrimitive).
 */
template <DerivedFromOrbital PrimitiveType>
class ContractedOrbital final : public Orbital {
  public:
    ContractedOrbital() = default;

    /**
     * @param size Number of primitives for the linear combination.
     */
    ContractedOrbital(size_t size) : m_coefficients(size), m_primitives(size) {}

    /**
     * @param coefficients Weights for each primitive in the linear combination.
     * @param primitives Primitives basis functions.
     * @throws std::length_error if the number of coefficients and primitives
     * are not equal.
     */
    ContractedOrbital(std::vector<double> &coefficients,
                      std::vector<PrimitiveType> &primitives);

    ContractedOrbital(const ContractedOrbital<PrimitiveType> &other);

    ContractedOrbital(const ContractedOrbital<PrimitiveType> &&other);

    /**
     * @brief Allocate memory for primitives and coefficients.
     *
     * @param size Number of primitives for the linear combination.
     */
    void reserve(std::size_t size);

    /**
     * @brief Add a primitive to the linear combination.
     *
     * @param coefficient Weight of the primitive in the linear combination.
     * @param primitive Primitive to add.
     */
    void add_primitive(double coefficient, const PrimitiveType &primitive);

    /**
     * @brief Add a primitive to the linear combination.
     *
     * @tparam Args
     * @param coefficient Weight of the primitive in the linear combination.
     * @param args Arguments to pass to the PrimitiveType constructor.
     */
    template <typename... Args>
    void add_primitive(double coefficient, Args &&...args);

    void update_normalization();

    inline size_t size() const;

    inline const PrimitiveType &get_primitive(int i) const;
    inline double get_coefficient(int i) const;

    void set_position(Eigen::Vector3d position) override {
        m_position = position;

        for (auto &primitive : m_primitives) {
            primitive.set_position(position);
        }
    }

  private:
    std::vector<double> m_coefficients;
    std::vector<PrimitiveType> m_primitives;
};

extern template class ContractedOrbital<GaussianPrimitive>;
extern template class ContractedOrbital<SlaterPrimitive>;

typedef ContractedOrbital<GaussianPrimitive> ContractedGaussian;
typedef ContractedOrbital<SlaterPrimitive> ContractedSlater;

typedef ContractedOrbital<GaussianPrimitive> CGTO;
typedef ContractedOrbital<SlaterPrimitive> CSTO;

template <DerivedFromOrbital PrimitiveType>
double overlap_integral(const ContractedOrbital<PrimitiveType> &,
                        const ContractedOrbital<PrimitiveType> &);

template <DerivedFromOrbital PrimitiveType>
double laplacian_integral(const ContractedOrbital<PrimitiveType> &,
                          const ContractedOrbital<PrimitiveType> &);

template <DerivedFromOrbital PrimitiveType>
double electron_nucleus_integral(const ContractedOrbital<PrimitiveType> &,
                                 const ContractedOrbital<PrimitiveType> &,
                                 const Eigen::Vector3d &);

template <DerivedFromOrbital PrimitiveType>
double electron_electron_integral(const ContractedOrbital<PrimitiveType> &,
                                  const ContractedOrbital<PrimitiveType> &,
                                  const ContractedOrbital<PrimitiveType> &,
                                  const ContractedOrbital<PrimitiveType> &);