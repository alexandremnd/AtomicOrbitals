#pragma once

#include <iostream>
#include <vector>

#include "Atom/atom.hpp"
#include "BasisSet/contracted_orbital.hpp"
#include "BasisSet/gaussian_primitive.hpp"
#include "BasisSet/orbital.hpp"
#include "BasisSet/slater_primitive.hpp"
#include "concepts.hpp"

/**
 * @brief Provides a class to model molecule that can be used in Hartree-Fock.
 * Provides many methods to easen iteration over bond length, orbital
 * optimizations and many more.
 *
 * @tparam OrbitalType The type of orbitals to use (SlaterPrimitive,
 * GaussianPrimitive, ...).
 */
template <DerivedFromOrbital OrbitalType> class Molecule {
  public:
    Molecule() = default;
    Molecule(const Molecule &molecule)
        : m_atoms(molecule.m_atoms),
          m_nucleus_charge(molecule.m_nucleus_charge),
          m_positions(molecule.m_positions), m_orbitals(molecule.m_orbitals) {};
    Molecule(Molecule &&molecule)
        : m_atoms(std::move(molecule.m_atoms)),
          m_nucleus_charge(std::move(molecule.m_nucleus_charge)),
          m_positions(std::move(molecule.m_positions)),
          m_orbitals(std::move(molecule.m_orbitals)) {};

    friend std::ostream &operator<<(std::ostream &os,
                                    const Molecule<OrbitalType> &molecule);

    /**
     * @brief Adds an atom to the current molecule.
     * @note The atom is copied into the molecule. Updating passed atom instance
     * will not affect the molecule.
     *
     * @param atom Atom to add
     */
    void add_atom(const Atom<OrbitalType> &atom);

    /**
     * @brief Adds an atom to the current molecule.
     * @note The atom is moved into the molecule, do not use the passed instance
     * as it is now in a null state.
     *
     * @param atom Atom to add.
     */
    void add_atom(Atom<OrbitalType> &&atom);

    /**
     * @brief Adds an atom to the current molecule
     *
     * @param nucleus_charge
     * @param position
     */
    void add_atom(int nucleus_charge,
                  const Eigen::Vector3d &position = Eigen::Vector3d::Zero());

    /**
     * @brief Adds an orbital to the molecule.
     *
     * @param orbital
     */
    void add_orbital(OrbitalType orbital) {
        m_orbitals.push_back(std::ref(orbital));
    }

    template <typename... Args>
    void add_orbital(Args &&...args)
        requires std::is_constructible_v<OrbitalType, Args...>
    {
        m_orbitals.emplace_back(std::forward<Args>(args)...);
    }

    OrbitalType &get_orbital(size_t atom, size_t orbital) const {
        return m_atoms[atom].get_orbital(orbital);
    }

    OrbitalType &get_orbital(size_t orbital) const {
        return m_orbitals[orbital];
    }
    Atom<OrbitalType> &get_atom(size_t atom) const { return m_atoms[atom]; }

  private:
    std::vector<Atom<OrbitalType>> m_atoms;
    std::vector<int> m_nucleus_charge;
    std::vector<Eigen::Vector3d> m_positions;
    std::vector<std::reference_wrapper<OrbitalType>> m_orbitals;

    // reference_wrapper is enough here as the lifetime of the orbitals are tied
    // to the atoms owned by m_atoms.
};

DECLARE_EXTERN_TEMPLATE(Molecule)