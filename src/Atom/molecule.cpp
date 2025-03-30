#include "Atom/molecule.hpp"

// ===============================
// System interface implementation
// ===============================
size_t Molecule::size() const { return m_orbitals.size(); }

double Molecule::overlap(size_t i, size_t j) const {
    return overlap_integral(m_orbitals[i].get(), m_orbitals[j].get());
}

double Molecule::kinetic(size_t i, size_t j) const {
    return -0.5 * laplacian_integral(m_orbitals[i].get(), m_orbitals[j].get());
}

double Molecule::electron_nucleus(size_t i, size_t j) const {
    double attraction = 0.0;

    for (size_t k = 0; k < m_atoms.size(); k++) {
        attraction +=
            -m_atoms[k]->Z() *
            electron_nucleus_integral(m_orbitals[i].get(), m_orbitals[j].get(),
                                      m_atoms[k]->position());
    }

    return attraction;
}

double Molecule::electron_electron(size_t i, size_t j, size_t k,
                                   size_t l) const {
    return electron_electron_integral(m_orbitals[i].get(), m_orbitals[j].get(),
                                      m_orbitals[k].get(), m_orbitals[l].get());
}

double Molecule::nucleus_repulsion() const {
    double repulsion = 0.0;

    for (size_t i = 0; i < m_atoms.size(); i++) {
        for (size_t j = i + 1; j < m_atoms.size(); j++) {
            repulsion += m_atoms[i]->Z() * m_atoms[j]->Z() / distance(i, j);
        }
    }
    return repulsion;
}

// ====================================
// Misceleaneous Functions
// ====================================

std::ostream &operator<<(std::ostream &os, const Molecule &mol) {
    for (const auto &atom : mol.m_atoms) {
        os << static_cast<Element>(atom->Z()) << " "
           << atom->position().transpose() << "\n";
    }
    os << "****\n";

    for (const auto &orb : mol.m_orbitals) {
        int i = 0;
        for (const auto &sub_orb : orb.get().get_primitives()) {
            os << orb.get().get_coefficient(i) << " " << sub_orb.constant();
            os << " " << sub_orb.position().transpose();
            os << " " << sub_orb.alpha() << " " << sub_orb.x_exponent() << " ";
            os << sub_orb.y_exponent() << " " << sub_orb.z_exponent() << "\n";
            i++;
        }
        os << "****\n";
    }

    return os;
}