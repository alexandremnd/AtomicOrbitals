#pragma once

#include <iostream>
#include <memory>
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
    Molecule(const Molecule &molecule) : m_atoms(molecule.m_atoms) {};
    Molecule(Molecule &&molecule) : m_atoms(std::move(molecule.m_atoms)) {}
    friend std::ostream &operator<<(std::ostream &os,
                                    const Molecule<OrbitalType> &molecule);

    void add_atom(std::shared_ptr<Atom<OrbitalType>> atom) { m_atoms.push_back(atom); }

    Atom<OrbitalType> &get_atom(size_t atom) const { return m_atoms[atom]; }

    const double distance(size_t i, size_t j) const {
        return (m_atoms[i].position() - m_atoms[j].position()).norm();
    }

  private:
    std::vector<std::shared_ptr<Atom<OrbitalType>>> m_atoms;
};


template <DerivedFromOrbital OrbitalType>
std::ostream &operator<<(std::ostream &os,
                         const Molecule<OrbitalType> &molecule) {
    os << "============= Molecule Configuration =============\n";
    os << "Number of atoms: " << molecule.m_atoms.size() << "\n";
    os << "Number of orbitals: " << molecule.m_orbitals.size() << "\n";
    os << "Atoms in molecule: " << "\n";
    for (size_t i = 0; i < molecule.m_atoms.size(); i++) {
        os << "\tAtom no." << i << " - Z = " << molecule.get_atom(i) << "\n";
        os << "\t\t- Position: " << molecule.get_atom(i).position() << "\n";
        os << "\t\t- Orbitals: " << "\n";
    }
    os << "==================================================\n";
}

DECLARE_EXTERN_TEMPLATE(Molecule)