#pragma once

#include "Atom/atom.hpp"
#include <vector>

template <typename T>
class Molecule {
    public:
        Molecule() {};
        Molecule(std::vector<Atom<T>>& atoms) : m_atoms(atoms) {};
        Molecule(const Molecule& molecule) : m_atoms(molecule.m_atoms) {};
        Molecule(const Molecule&& molecule) : m_atoms(std::move(molecule.m_atoms)) {};

        void add_atom(const Atom<T>& atom) {
            m_atoms.emplace_back(atom);
        }

        inline const int n_atoms() const { return m_atoms.size(); }
        inline const std::vector<Atom<T>>& atoms() const { return m_atoms; }
        inline const Atom<T>& atom(size_t i) const { return m_atoms[i]; }
        inline const Atom<T>& operator()(size_t i) const { return m_atoms[i]; }

    private:
        std::vector<Atom<T>> m_atoms;
};