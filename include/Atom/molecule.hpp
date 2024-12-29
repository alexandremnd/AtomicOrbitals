#pragma once

#include "Atom/atom.hpp"
#include <memory>
#include <vector>

template <typename T>
class Molecule {
    public:
        Molecule() {};
        Molecule(std::vector<Atom<T>>& atoms) : m_atoms(atoms) {};
        Molecule(const Molecule& molecule) : m_atoms(molecule.m_atoms) {};
        Molecule(const Molecule&& molecule) : m_atoms(std::move(molecule.m_atoms)) {};

        /**
         * Provides a copyless way to add an atom to the molecule (useful for many atoms molecules)
         */
        void add_atom(std::unique_ptr<Atom<T>>&& atom) {
            m_atoms.push_back(std::move(atom));
        }

        void add_atom(const Atom<T>& atom) {
            m_atoms.push_back(std::make_unique<Atom<T>>(atom));
        }

        inline const int n_atoms() const { return m_atoms.size(); }
        inline const std::vector<std::unique_ptr<Atom<T>>>& atoms() const { return m_atoms; }
        inline const std::unique_ptr<Atom<T>>& atom(size_t i) const { return m_atoms[i]; }
        inline const std::unique_ptr<Atom<T>>& operator()(size_t i) const { return m_atoms[i]; }

    private:
        std::vector<std::unique_ptr<Atom<T>>> m_atoms;
};