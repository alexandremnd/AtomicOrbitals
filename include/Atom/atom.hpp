#pragma once

#include <stdexcept>
#include <type_traits>
#include <vector>

#include "Atom/system.hpp"
#include "Orbitals/contracted_gaussian.hpp"
#include "Atom/atom_list.hpp"

#include "Eigen/Dense"

/**
 * @brief Class representing an atom (if not obvious ...).
 *
 * This is the privileged class to build any system that may be forwarded to
 * Hartree-Fock for computation. Molecules are built upon atoms.
 */
class Atom : public System {
  public:
    /**
     * @brief Builds an atom at the given position.
     * @param Z Number of protons in the nucleus.
     * @param position Position of the nucleus.
     * @throw std::invalid_argument if Z < 1.
     */
    Atom(int Z, Eigen::Vector3d position = {0, 0, 0})
        : m_Z(Z), m_position(position) {
        if (Z < 1) {
            throw std::invalid_argument(
                "Atom: The atomic number must be greater than 0.");
        }
    };

    /**
     * @brief Builds an atom at the given position and loads corresponding
     * basis.
     * @note Basis name is case sensitive, by convention, all basis file are in
     * lowercase.
     *
     * @param elt Element to build (He, Li, ...)
     * @param basis_name Basis type to use (sto-6g, sto, ugbs, ...)
     * @param position Position of the atom
     */
    Atom(Element elt, std::string basis_name,
         Eigen::Vector3d position = {0, 0, 0});

    Atom(Atom &atom)
        : m_Z(atom.Z()), m_position(atom.m_position),
          m_orbitals(atom.m_orbitals) {};
    Atom(Atom &&atom)
        : m_Z(atom.Z()), m_position(atom.m_position),
          m_orbitals(std::move(atom.m_orbitals)) {};

    void print_info() const;

    /**
     * @param orbital The orbital to add to the atom.
     * @note The orbital will be copied in this case.
     */
    void add_orbital(const ContractedGaussian &orbital);

    /**
     * @param orbital The orbital to add to the atom.
     * @note The orbital will be moved in this case. Passed orbital will be
     * invalid after this call.
     */
    void add_orbital(ContractedGaussian &&orbital);

    /**
     * @brief Constructs an orbital with T constructor
     * @param args Arguments to pass to the constructor of the orbital.
     */
    template <typename... Args>
    void add_orbital(Args &&...args)
        requires std::is_constructible_v<ContractedGaussian, Args...>
    {
        m_orbitals.emplace_back(std::forward<Args>(args)...);
    }

    /**
     * @brief Adds a contracted gaussian orbital (s type) to the atom.
     *
     * The contracted gaussian orbital will be a linear combination of
     * weight.size() primitive gaussians.
     *
     * @throw std::invalid_argument if (weight.size() != decay.size())
     *
     * @param weight Coefficient of the primitive gaussians in the linear
     * combination.
     * @param decay Exponential decay rate of the primitive gaussians.
     */
    void add_gaussian_orbital_stype(const std::vector<double> &weight,
                                    const std::vector<double> &decay);

    /**
     * @brief Adds a contracted gaussian orbital (p type) to the atom.
     *
     * The contracted gaussian orbital will be a linear combination of
     * weight.size() primitive gaussians.
     *
     * @throw std::invalid_argument if (weight.size() != decay.size())
     * @param weight Coefficient of the primitive gaussians in the linear
     * combination.
     * @param decay Exponential decay rate of the primitive gaussians.
     */
    void add_gaussian_orbital_ptype(const std::vector<double> &weight,
                                    const std::vector<double> &decay);

    /**
     * @brief Adds a contracted gaussian orbital (d type) to the atom.
     *
     * The contracted gaussian orbital will be a linear combination of
     * weight.size() primitive gaussians.
     *
     * @throw std::invalid_argument if (weight.size() != decay.size())
     * @param weight Coefficient of the primitive gaussians in the linear
     * combination.
     * @param decay Exponential decay rate of the primitive gaussians.
     */
    void add_gaussian_orbital_dtype(const std::vector<double> &weight,
                                    const std::vector<double> &decay);

    /**
     * @brief Adds a contracted gaussian orbital (f type ie l = 3) to the atom.
     *
     * The contracted gaussian orbital will be a linear combination of
     * weight.size() primitive gaussians.
     *
     * @throw std::invalid_argument if (weight.size() != decay.size())
     * @param weight Coefficient of the primitive gaussians in the linear
     * combination.
     * @param decay Exponential decay rate of the primitive gaussians.
     */
    void add_gaussian_orbital_ftype(const std::vector<double> &weight,
                                    const std::vector<double> &decay);

    /**
     * @brief Sets the position of the nucleus.
     * @note This will also update the center of the orbitals if required.
     *
     * @param position New position of the nucleus.
     */
    void set_position(const Eigen::Vector3d &position);

    inline int Z() { return m_Z; }
    inline Eigen::Vector3d position() { return m_position; }

    inline std::vector<ContractedGaussian> &get_orbitals() {
        return m_orbitals;
    }
    inline ContractedGaussian &get_orbital(size_t i) { return m_orbitals[i]; }

    // ===============================
    // System interface implementation
    // ===============================
    size_t size() const override;
    double overlap(size_t i, size_t j) const override;

    double kinetic(size_t i, size_t j) const override;

    double electron_nucleus(size_t i, size_t j) const override;

    double electron_electron(size_t i, size_t j, size_t k,
                             size_t l) const override;
    double nucleus_repulsion() const override;

    friend std::ostream &operator<<(std::ostream &os, const Atom &orbital);

  private:
    int m_Z;
    Eigen::Vector3d m_position;
    std::vector<ContractedGaussian> m_orbitals;
};