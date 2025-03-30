#pragma once

#include <functional>
#include <memory>
#include <vector>
#include <iostream>

#include "Atom/atom.hpp"
#include "Atom/system.hpp"
#include "Orbitals/contracted_gaussian.hpp"

/**
 * @brief Provides a class to model molecule that can be used in Hartree-Fock.
 *
 */
class Molecule : public System {
  public:
    Molecule() = default;
    Molecule(const Molecule &molecule) : m_atoms(molecule.m_atoms) {};
    Molecule(Molecule &&molecule) : m_atoms(std::move(molecule.m_atoms)) {}

    void add_atom(std::shared_ptr<Atom> atom) {
        m_atoms.push_back(atom);

        for (ContractedGaussian &orbital : atom->get_orbitals()) {
            m_orbitals.emplace_back(orbital);
        }
    }

    Atom &get_atom(size_t i) { return *m_atoms[i]; }
    ContractedGaussian &get_orbital(size_t i) { return m_orbitals[i]; }
    std::vector<std::reference_wrapper<ContractedGaussian>> &get_orbitals() {
        return m_orbitals;
    }

    const double distance(size_t i, size_t j) const {
        return (m_atoms[i]->position() - m_atoms[j]->position()).norm();
    }

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

    friend std::ostream &operator<<(std::ostream &os, const Molecule &atom);

  private:
    std::vector<std::shared_ptr<Atom>> m_atoms;
    std::vector<std::reference_wrapper<ContractedGaussian>> m_orbitals;
};
