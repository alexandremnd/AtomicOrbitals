#pragma once

#include <functional>
#include <iostream>
#include <memory>
#include <vector>

#include "Atom/atom.hpp"
#include "Atom/system.hpp"
#include "Orbitals/contracted_orbital.hpp"
#include "Orbitals/gaussian_primitive.hpp"
#include "Orbitals/slater_primitive.hpp"
#include "concepts.hpp"

/**
 * @brief Provides a class to model molecule that can be used in Hartree-Fock.
 * Provides many methods to easen iteration over bond length, orbital
 * optimizations and many more.
 *
 * @tparam OrbitalType The type of orbitals to use (SlaterPrimitive,
 * GaussianPrimitive, ...).
 */
template <DerivedFromOrbital OrbitalType> class Molecule : public System {
  public:
    Molecule() = default;
    Molecule(const Molecule &molecule) : m_atoms(molecule.m_atoms){};
    Molecule(Molecule &&molecule) : m_atoms(std::move(molecule.m_atoms)) {}

    void add_atom(std::shared_ptr<Atom<OrbitalType>> atom) {
        m_atoms.push_back(atom);

        for (OrbitalType &orbital : atom->get_orbitals()) {
            m_orbitals.emplace_back(orbital);
        }
    }

    Atom<OrbitalType> &get_atom(size_t i) { return *m_atoms[i]; }
    OrbitalType &get_orbital(size_t i) { return m_orbitals[i]; }

    const double distance(size_t i, size_t j) const {
        return (m_atoms[i]->position() - m_atoms[j]->position()).norm();
    }

    // System interface implementation
    size_t size() const override { return m_orbitals.size(); }

    double overlap(size_t i, size_t j) const override {
        return overlap_integral(m_orbitals[i].get(), m_orbitals[j].get());
    }

    double kinetic(size_t i, size_t j) const override {
        return -0.5 *
               laplacian_integral(m_orbitals[i].get(), m_orbitals[j].get());
    }

    double electron_nucleus(size_t i, size_t j) const override {
        double attraction = 0.0;

        for (size_t k = 0; k < m_atoms.size(); k++) {
            attraction += -m_atoms[k]->Z() *
                          electron_nucleus_integral(m_orbitals[i].get(),
                                                    m_orbitals[j].get(),
                                                    m_atoms[k]->position());
        }

        return attraction;
    }

    double electron_electron(size_t i, size_t j, size_t k,
                             size_t l) const override {
        return electron_electron_integral(
            m_orbitals[i].get(), m_orbitals[j].get(), m_orbitals[k].get(),
            m_orbitals[l].get());
    }

    double nuclear_energy() const override {
        double repulsion = 0.0;

        for (size_t i = 0; i < m_atoms.size(); i++) {
            for (size_t j = i + 1; j < m_atoms.size(); j++) {
                repulsion += m_atoms[i]->Z() * m_atoms[j]->Z() / distance(i, j);
            }
        }

        return repulsion;
    }

  private:
    std::vector<std::shared_ptr<Atom<OrbitalType>>> m_atoms;
    std::vector<std::reference_wrapper<OrbitalType>> m_orbitals;
};

DECLARE_EXTERN_TEMPLATE(Molecule)