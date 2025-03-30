#include <iostream>
#include <stdexcept>

#include "Atom/atom.hpp"
#include "Orbitals/contracted_gaussian.hpp"
#include "basis_parser.hpp"

Atom::Atom(Element elt, std::string basis_name, Eigen::Vector3d position) {
    m_position = position;
    m_Z = static_cast<int>(elt);
    parse_basis(elt, basis_name, *this);
}

void Atom::print_info() const {
    std::cout << "============= Atom Configuration =============" << std::endl;
    std::cout << "Atomic number: " << m_Z << std::endl;
    std::cout << "Position: " << m_position.transpose() << std::endl;
    std::cout << "Orbitals count: " << m_orbitals.size() << std::endl;
    std::cout << "==============================================" << std::endl;
}

void Atom::add_orbital(const ContractedGaussian &orbital) {
    m_orbitals.push_back(orbital);
}

void Atom::add_orbital(ContractedGaussian &&orbital) {
    m_orbitals.push_back(std::move(orbital));
}

// ====================================
// Contracted Gaussian Orbital Functions
// ====================================

void Atom::add_gaussian_orbital_stype(const std::vector<double> &weight,
                                      const std::vector<double> &decay) {
    if (weight.size() != decay.size()) {
        throw std::invalid_argument(
            "Atom: The weight and decay vectors must have the same size.");
    }

    auto cg = ContractedGaussian();
    for (size_t i = 0; i < weight.size(); i++) {
        cg.add_primitive(weight[i], 0, 0, 0, decay[i], m_position);
    }

    m_orbitals.push_back(std::move(cg));
}

void Atom::add_gaussian_orbital_ptype(const std::vector<double> &weight,
                                      const std::vector<double> &decay) {
    if (weight.size() != decay.size()) {
        throw std::invalid_argument(
            "Atom: The weight and decay vectors must have the same size.");
    }

    auto cg_x = ContractedGaussian();
    auto cg_y = ContractedGaussian();
    auto cg_z = ContractedGaussian();

    for (size_t i = 0; i < weight.size(); i++) {
        cg_x.add_primitive(weight[i], 1, 0, 0, decay[i], m_position);
        cg_y.add_primitive(weight[i], 0, 1, 0, decay[i], m_position);
        cg_z.add_primitive(weight[i], 0, 0, 1, decay[i], m_position);
    }

    m_orbitals.push_back(std::move(cg_x));
    m_orbitals.push_back(std::move(cg_y));
    m_orbitals.push_back(std::move(cg_z));
}

void Atom::add_gaussian_orbital_dtype(const std::vector<double> &weight,
                                      const std::vector<double> &decay) {
    if (weight.size() != decay.size()) {
        throw std::invalid_argument(
            "Atom: The weight and decay vectors must have the same size.");
    }

    auto cg_xx = ContractedGaussian();
    auto cg_yy = ContractedGaussian();
    auto cg_zz = ContractedGaussian();
    auto cg_xy = ContractedGaussian();
    auto cg_xz = ContractedGaussian();
    auto cg_yz = ContractedGaussian();

    for (size_t i = 0; i < weight.size(); i++) {
        cg_xx.add_primitive(weight[i], 2, 0, 0, decay[i], m_position);
        cg_yy.add_primitive(weight[i], 0, 2, 0, decay[i], m_position);
        cg_zz.add_primitive(weight[i], 0, 0, 2, decay[i], m_position);
        cg_xy.add_primitive(weight[i], 1, 1, 0, decay[i], m_position);
        cg_xz.add_primitive(weight[i], 1, 0, 1, decay[i], m_position);
        cg_yz.add_primitive(weight[i], 0, 1, 1, decay[i], m_position);
    }

    m_orbitals.push_back(std::move(cg_xx));
    m_orbitals.push_back(std::move(cg_yy));
    m_orbitals.push_back(std::move(cg_zz));
    m_orbitals.push_back(std::move(cg_xy));
    m_orbitals.push_back(std::move(cg_xz));
    m_orbitals.push_back(std::move(cg_yz));
}

void Atom::add_gaussian_orbital_ftype(const std::vector<double> &weight,
                                      const std::vector<double> &decay) {
    if (weight.size() != decay.size()) {
        throw std::invalid_argument(
            "Atom: The weight and decay vectors must have the same size.");
    }

    auto cg_xxx = ContractedGaussian();
    auto cg_yyy = ContractedGaussian();
    auto cg_zzz = ContractedGaussian();
    auto cg_xxy = ContractedGaussian();
    auto cg_xxz = ContractedGaussian();
    auto cg_xyy = ContractedGaussian();
    auto cg_yyz = ContractedGaussian();
    auto cg_xzz = ContractedGaussian();
    auto cg_yzz = ContractedGaussian();
    auto cg_xyz = ContractedGaussian();

    for (size_t i = 0; i < weight.size(); i++) {
        cg_xxx.add_primitive(weight[i], 3, 0, 0, decay[i], m_position);
        cg_yyy.add_primitive(weight[i], 0, 3, 0, decay[i], m_position);
        cg_zzz.add_primitive(weight[i], 0, 0, 3, decay[i], m_position);
        cg_xxy.add_primitive(weight[i], 2, 1, 0, decay[i], m_position);
        cg_xxz.add_primitive(weight[i], 2, 0, 1, decay[i], m_position);
        cg_xyy.add_primitive(weight[i], 1, 2, 0, decay[i], m_position);
        cg_yyz.add_primitive(weight[i], 0, 2, 1, decay[i], m_position);
        cg_xzz.add_primitive(weight[i], 1, 0, 2, decay[i], m_position);
        cg_yzz.add_primitive(weight[i], 0, 1, 2, decay[i], m_position);
        cg_xyz.add_primitive(weight[i], 1, 1, 1, decay[i], m_position);
    }

    m_orbitals.push_back(std::move(cg_xxx));
    m_orbitals.push_back(std::move(cg_yyy));
    m_orbitals.push_back(std::move(cg_zzz));
    m_orbitals.push_back(std::move(cg_xxy));
    m_orbitals.push_back(std::move(cg_xxz));
    m_orbitals.push_back(std::move(cg_xyy));
    m_orbitals.push_back(std::move(cg_yyz));
    m_orbitals.push_back(std::move(cg_xzz));
    m_orbitals.push_back(std::move(cg_yzz));
    m_orbitals.push_back(std::move(cg_xyz));
}

// ===============================
// System interface implementation
// ===============================

size_t Atom::size() const { return m_orbitals.size(); }

double Atom::overlap(size_t i, size_t j) const {
    return overlap_integral(m_orbitals[i], m_orbitals[j]) *
           m_orbitals[i].constant() * m_orbitals[j].constant();
}

double Atom::kinetic(size_t i, size_t j) const {
    return -0.5 * laplacian_integral(m_orbitals[i], m_orbitals[j]) *
           m_orbitals[i].constant() * m_orbitals[j].constant();
}

double Atom::electron_nucleus(size_t i, size_t j) const {
    return -m_Z *
           electron_nucleus_integral(m_orbitals[i], m_orbitals[j], m_position) *
           m_orbitals[i].constant() * m_orbitals[j].constant();
}

double Atom::electron_electron(size_t i, size_t j, size_t k, size_t l) const {
    return electron_electron_integral(m_orbitals[i], m_orbitals[j],
                                      m_orbitals[k], m_orbitals[l]) *
           m_orbitals[i].constant() * m_orbitals[j].constant() *
           m_orbitals[k].constant() * m_orbitals[l].constant();
}

double Atom::nucleus_repulsion() const { return 0.; }

// ====================================
// Misceleaneous Functions
// ====================================

std::ostream &operator<<(std::ostream &os, const Atom &atom) {
    os << static_cast<Element>(atom.m_Z) << " " << atom.m_position.transpose()
       << "\n";
    os << "****\n";

    for (const auto &orb : atom.m_orbitals) {
        int i = 0;
        for (const auto &sub_orb : orb.get_primitives()) {
            os << orb.get_coefficient(i) << " " << sub_orb.constant();
            os << " " << sub_orb.position().transpose();
            os << " " << sub_orb.alpha() << " " << sub_orb.x_exponent() << " ";
            os << sub_orb.y_exponent() << " " << sub_orb.z_exponent() << "\n";
            i++;
        }
        os << "****\n";
    }

    return os;
}

void Atom::set_position(const Eigen::Vector3d &position) {
    m_position = position;

    for (ContractedGaussian &orbital : m_orbitals) {
        orbital.set_position(m_position);
    }
}