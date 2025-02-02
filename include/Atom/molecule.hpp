#include "Atom/molecule.interface.hpp"
#include "concepts.hpp"

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

template <DerivedFromOrbital OrbitalType>
void Molecule<OrbitalType>::add_atom(const Atom<OrbitalType> &atom) {
    m_atoms.push_back(atom);
    m_nucleus_charge.push_back(atom.Z());
    m_positions.push_back(atom.position());

    for (const auto &orbital : atom.get_orbitals()) {
        m_orbitals.push_back(std::ref(orbital));
    }
}

template <DerivedFromOrbital OrbitalType>
void Molecule<OrbitalType>::add_atom(Atom<OrbitalType> &&atom) {
    m_nucleus_charge.push_back(atom.Z());
    m_positions.push_back(atom.position());

    for (const auto &orbital : atom.get_orbitals()) {
        m_orbitals.push_back(std::ref(orbital));
    }

    m_atoms.push_back(std::move(atom));
}

template <DerivedFromOrbital OrbitalType>
void Molecule<OrbitalType>::add_atom(int nucleus_charge,
                                     const Eigen::Vector3d &position) {
    m_nucleus_charge.push_back(nucleus_charge);
    m_positions.push_back(position);
    m_atoms.emplace_back(nucleus_charge, position);
}