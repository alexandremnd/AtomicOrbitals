#pragma once

#include <iostream>

#include "Atom/atom.interface.hpp" // IWYU pragma: export
#include "Orbitals/contracted_orbital.hpp"
#include "parse_basis.hpp"
#include "concepts.hpp"

template <DerivedFromOrbital OrbitalType>
Atom<OrbitalType>::Atom(Element elt, std::string basis_name,
                        Eigen::Vector3d position) {
    m_Z = elt;
    m_position = position;
    parse_basis(elt, basis_name, *this);
}

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::print_info() const {
    std::cout << "============= Atom Configuration =============" << std::endl;
    std::cout << "Atomic number: " << m_Z << std::endl;
    std::cout << "Position: " << m_position.transpose() << std::endl;
    std::cout << "Orbitals count: " << m_orbitals.size() << std::endl;
    std::cout << "==============================================" << std::endl;
}

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::add_orbital(const OrbitalType &orbital) {
    m_orbitals.push_back(orbital);
}

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::add_orbital(OrbitalType &&orbital) {
    m_orbitals.push_back(std::move(orbital));
}

template <DerivedFromOrbital OrbitalType>
template <typename... Args>
void Atom<OrbitalType>::add_orbital(Args &&...args)
    requires std::is_constructible_v<OrbitalType, Args...>
{
    m_orbitals.emplace_back(std::forward<Args>(args)...);
}

// ====================================
// Contracted Slater Orbital Functions
// ====================================

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::add_slater_orbital(const std::vector<double> &weight,
                                           const std::vector<double> &n,
                                           const std::vector<double> &l,
                                           const std::vector<double> &m,
                                           const std::vector<double> &decay)
    requires std::is_same_v<OrbitalType, ContractedSlater>
{
    if (weight.size() != decay.size()) {
        throw std::invalid_argument(
            "Atom: The weight and decay vectors must have the same size.");
    }

    auto cs = ContractedSlater();
    cs.reserve(weight.size());

    for (size_t i = 0; i < weight.size(); i++) {
        cs.add_primitive(weight[i], n[i], l[i], m[i], decay[i]);
    }

    m_orbitals.push_back(std::move(cs));
}

// ====================================
// Contracted Gaussian Orbital Functions
// ====================================

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::add_gaussian_orbital_stype(
    const std::vector<double> &weight, const std::vector<double> &decay)
    requires std::is_same_v<OrbitalType, ContractedGaussian>
{
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

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::add_gaussian_orbital_ptype(
    const std::vector<double> &weight, const std::vector<double> &decay)
    requires std::is_same_v<OrbitalType, ContractedGaussian>
{
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

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::add_gaussian_orbital_dtype(
    const std::vector<double> &weight, const std::vector<double> &decay)
    requires std::is_same_v<OrbitalType, ContractedGaussian>
{
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

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::add_gaussian_orbital_ftype(
    const std::vector<double> &weight, const std::vector<double> &decay)
    requires std::is_same_v<OrbitalType, ContractedGaussian>
{
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

// ====================================
// Misceleaneous Functions
// ====================================

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::set_position(const Eigen::Vector3d &position) {
    m_position = position;

    for (OrbitalType &orbital : m_orbitals) {
        orbital.set_position(m_position);
    }
}