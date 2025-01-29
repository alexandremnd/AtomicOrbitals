#pragma once

#include <memory>
#include <type_traits>
#include <vector>
#include <stdexcept>

#include "BasisSet/orbital.hpp"
#include "BasisSet/slater_primitive.hpp"
#include "BasisSet/slater_contracted.hpp"
#include "BasisSet/gaussian_contracted.hpp"

#include "concepts.hpp"
#include "Eigen/Dense"

template <DerivedFromOrbital OrbitalType>
class Atom {
public:
    /**
    * Throws an invalid_argument exception if Z < 1.
    *
    * @param Z Number of protons in the nucleus.
    * @param position Position of the nucleus.
    */
    Atom(int Z, Eigen::Vector3d position) : m_Z(Z), m_position(position) {
        if (Z < 1) {
            throw std::invalid_argument("Atom: The atomic number must be greater than 0.");
        }
    };

    /**
    * @brief Builds an atom at the origin.
    *
    * Throws an invalid_argument exception if Z < 1.
    *
    * @param Z Number of protons in the nucleus.
    */
    Atom(int Z) : Atom(Z, Eigen::Vector3d{0., 0., 0.}) {};

    /**
    * Throws an invalid_argument exception if Z < 1.
    *
    * @param Z Number of protons in the nucleus.
    * @note Orbitals are not copied. Only pointers are copied.
    */
    Atom(int Z, float x, float y, float z) : Atom(Z, Eigen::Vector3d{x, y, z}) {};
    Atom(Atom& atom) : m_Z(atom.Z()), m_position(atom.m_position), m_orbitals(atom.m_orbitals) {};
    Atom(Atom&& atom) : m_Z(atom.Z()), m_position(atom.m_position), m_orbitals(std::move(atom.m_orbitals)) {};

    void print_info() const;

    /**
    * @param orbital The orbital to add to the atom.
    * @note The orbital will be copied in this case.
    */
    void add_orbital(const OrbitalType& orbital);

    /**
    * @param orbital The orbital to add to the atom.
    * @note The orbital will be moved in this case. The instance will own the orbital.
    */
    void add_orbital(OrbitalType&& orbital);

    /**
     * @brief Constructs an orbital with T constructor
     * @param args
     */
    template <typename... Args>
    void add_orbital(Args&&... args) requires std::is_constructible_v<OrbitalType, Args...> {
        m_orbitals->emplace_back(std::forward<Args>(args)...);
    }

    /**
    * @param n The principal quantum number. (0 < n)
    * @param l The azimuthal quantum number. (0 <= l < n)
    * @param m The magnetic quantum number. (-l <= m <= l)
    * @param alpha The radial exponential decay rate. (alpha > 0)
    */
    void add_slater_orbital(const int n, const int l, const int m, const double alpha) requires std::is_same_v<OrbitalType, SlaterPrimitive>;

    /**
    * @param weight The principal quantum number. (0 < n)
    * @param n The principal quantum number. (0 < n)
    * @param l The azimuthal quantum number. (0 <= l < n)
    * @param m The magnetic quantum number. (-l <= m <= l)
    * @param decay The exponential decay rate of the primitive gaussians (alpha > 0)
    */
    void add_contracted_slater(const std::vector<double> &weight,
                               const std::vector<double> &n,
                               const std::vector<double> &l,
                               const std::vector<double> &m,
                               const std::vector<double> &decay) requires std::is_same_v<OrbitalType, ContractedSlater>;

    /**
    * @brief Adds a contracted gaussian orbital (s type) to the atom.
    *
    * The contracted gaussian orbital will be a linear combination of weight.size() primitive gaussians.
    *
    * Throws an invalid_argument exception if (weight.size() != decay.size())
    *
    * @param weight Coefficient of the primitive gaussians in the linear combination.
    * @param decay Exponential decay rate of the primitive gaussians.
    */
    void add_gaussian_orbital_stype(const std::vector<double>& weight, const std::vector<double>& decay) requires std::is_same_v<OrbitalType, ContractedGaussian>;

    /**
    * @brief Add a contracted gaussian orbital (p type) to the atom.
    *
    * The contracted gaussian orbital will be a linear combination of weight.size() primitive gaussians.
    *
    * Throws an invalid_argument exception if (weight.size() != decay.size())
    * @param weight Coefficient of the primitive gaussians in the linear combination.
    * @param decay Exponential decay rate of the primitive gaussians.
    */
    void add_gaussian_orbital_ptype(const std::vector<double>& weight, const std::vector<double>& decay) requires std::is_same_v<OrbitalType, ContractedGaussian>;

    /**
    * @brief Add a contracted gaussian orbital (d type) to the atom.
    *
    * The contracted gaussian orbital will be a linear combination of weight.size() primitive gaussians.
    *
    * Throws an invalid_argument exception if (weight.size() != decay.size())
    * @param weight Coefficient of the primitive gaussians in the linear combination.
    * @param decay Exponential decay rate of the primitive gaussians.
    */
    void add_gaussian_orbital_dtype(const std::vector<double>& weight, const std::vector<double>& decay) requires std::is_same_v<OrbitalType, ContractedGaussian>;

    /**
    * @brief Add a contracted gaussian orbital (f type ie l = 3) to the atom.
    *
    * The contracted gaussian orbital will be a linear combination of weight.size() primitive gaussians.
    *
    * Throws an invalid_argument exception if (weight.size() != decay.size())
    * @param weight Coefficient of the primitive gaussians in the linear combination.
    * @param decay Exponential decay rate of the primitive gaussians.
    */
    void add_gaussian_orbital_ftype(const std::vector<double>& weight, const std::vector<double>& decay) requires std::is_same_v<OrbitalType, ContractedGaussian>;

    inline int Z() { return m_Z; }
    inline Eigen::Vector3d position() { return m_position; }

    inline int orbitals_count() { return m_orbitals.size(); }
    inline const std::vector<std::shared_ptr<Orbital>>& get_orbitals() const { return m_orbitals; }
    inline const Orbital& get_orbital(size_t i) const { return *m_orbitals[i]; }
    inline const Orbital& operator()(size_t i) const { return *m_orbitals[i]; }

private:
    int m_Z;
    Eigen::Vector3d m_position;
    // std::vector<std::shared_ptr<OrbitalType>> m_orbitals;
    std::unique_ptr<std::vector<OrbitalType>> m_orbitals;
};

DECLARE_EXTERN_TEMPLATE(Atom)