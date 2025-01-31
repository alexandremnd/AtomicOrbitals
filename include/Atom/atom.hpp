#pragma once

#include <memory>
#include <iostream>

#include "Atom/atom.interface.hpp"
#include "BasisSet/gaussian_contracted.hpp"
#include "BasisSet/slater_primitive.hpp"
#include "BasisSet/slater_contracted.hpp"
#include "concepts.hpp"

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::print_info() const {
    std::cout << "============= Atom Configuration =============" << std::endl;
    std::cout << "Atomic number: " << m_Z << std::endl;
    std::cout << "Position: " << m_position.transpose() << std::endl;
    std::cout << "Orbitals count: " << m_orbitals.size() << std::endl;
    std::cout << "==============================================" << std::endl;
}

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::add_orbital(const OrbitalType& orbital) {
    m_orbitals.push_back(std::make_shared<OrbitalType>(orbital));
}

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::add_slater_orbital(int n, int l, int m, double alpha) requires std::is_same_v<OrbitalType, SlaterPrimitive>{
    m_orbitals.push_back(std::make_shared<SlaterPrimitive>(n, l, m, alpha));
}

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::add_contracted_slater(const std::vector<double>& weight,
                                                const std::vector<double> &n,
                                                const std::vector<double> &l,
                                                const std::vector<double> &m,
                                                const std::vector<double>& decay) requires std::is_same_v<OrbitalType, ContractedSlater>
{
    if (weight.size() != decay.size()) {
        throw std::invalid_argument("Atom: The weight and decay vectors must have the same size.");
    }

    auto cs = std::make_shared<ContractedSlater>();
    cs->reserve(weight.size());

    for (size_t i = 0; i < weight.size(); i++) {
        cs->add_primitive(weight[i], n[i], l[i], m[i], decay[i]);
    }

    m_orbitals.push_back(std::move(cs));
}

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::add_gaussian_orbital_stype(const std::vector<double>& weight, const std::vector<double>& decay) requires std::is_same_v<OrbitalType, ContractedGaussian> {
    if (weight.size() != decay.size()) {
        throw std::invalid_argument("Atom: The weight and decay vectors must have the same size.");
    }

    auto cg = std::make_shared<ContractedGaussian>();
    for (size_t i = 0; i < weight.size(); i++) {
        cg->add_primitive(weight[i], decay[i], 0, 0, 0, m_position);
    }

    m_orbitals.push_back(std::move(cg));
}

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::add_gaussian_orbital_ptype(const std::vector<double>& weight, const std::vector<double>& decay) requires std::is_same_v<OrbitalType, ContractedGaussian> {
    if (weight.size() != decay.size()) {
        throw std::invalid_argument("Atom: The weight and decay vectors must have the same size.");
    }

    auto cg_x = std::make_shared<ContractedGaussian>();
    auto cg_y = std::make_shared<ContractedGaussian>();
    auto cg_z = std::make_shared<ContractedGaussian>();

    for (size_t i = 0; i < weight.size(); i++) {
        cg_x->add_primitive(weight[i], decay[i], 1, 0, 0, m_position);
        cg_y->add_primitive(weight[i], decay[i], 0, 1, 0, m_position);
        cg_z->add_primitive(weight[i], decay[i], 0, 0, 1, m_position);
    }

    m_orbitals.push_back(std::move(cg_x));
    m_orbitals.push_back(std::move(cg_y));
    m_orbitals.push_back(std::move(cg_z));
}

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::add_gaussian_orbital_dtype(const std::vector<double>& weight, const std::vector<double>& decay) requires std::is_same_v<OrbitalType, ContractedGaussian> {
    if (weight.size() != decay.size()) {
        throw std::invalid_argument("Atom: The weight and decay vectors must have the same size.");
    }

    auto cg_xx = std::make_shared<ContractedGaussian>();
    auto cg_yy = std::make_shared<ContractedGaussian>();
    auto cg_zz = std::make_shared<ContractedGaussian>();
    auto cg_xy = std::make_shared<ContractedGaussian>();
    auto cg_xz = std::make_shared<ContractedGaussian>();
    auto cg_yz = std::make_shared<ContractedGaussian>();

    for (size_t i = 0; i < weight.size(); i++) {
        cg_xx->add_primitive(weight[i], decay[i], 2, 0, 0, m_position);
        cg_yy->add_primitive(weight[i], decay[i], 0, 2, 0, m_position);
        cg_zz->add_primitive(weight[i], decay[i], 0, 0, 2, m_position);
        cg_xy->add_primitive(weight[i], decay[i], 1, 1, 0, m_position);
        cg_xz->add_primitive(weight[i], decay[i], 1, 0, 1, m_position);
        cg_yz->add_primitive(weight[i], decay[i], 0, 1, 1, m_position);
    }

    m_orbitals.push_back(std::move(cg_xx));
    m_orbitals.push_back(std::move(cg_yy));
    m_orbitals.push_back(std::move(cg_zz));
    m_orbitals.push_back(std::move(cg_xy));
    m_orbitals.push_back(std::move(cg_xz));
    m_orbitals.push_back(std::move(cg_yz));
}

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::add_gaussian_orbital_ftype(const std::vector<double>& weight, const std::vector<double>& decay) requires std::is_same_v<OrbitalType, ContractedGaussian> {
    if (weight.size() != decay.size()) {
        throw std::invalid_argument("Atom: The weight and decay vectors must have the same size.");
    }

    auto cg_xxx = std::make_shared<ContractedGaussian>();
    auto cg_yyy = std::make_shared<ContractedGaussian>();
    auto cg_zzz = std::make_shared<ContractedGaussian>();
    auto cg_xxy = std::make_shared<ContractedGaussian>();
    auto cg_xxz = std::make_shared<ContractedGaussian>();
    auto cg_xyy = std::make_shared<ContractedGaussian>();
    auto cg_yyz = std::make_shared<ContractedGaussian>();
    auto cg_xzz = std::make_shared<ContractedGaussian>();
    auto cg_yzz = std::make_shared<ContractedGaussian>();
    auto cg_xyz = std::make_shared<ContractedGaussian>();

    for (size_t i = 0; i < weight.size(); i++) {
        cg_xxx->add_primitive(weight[i], decay[i], 3, 0, 0, m_position);
        cg_yyy->add_primitive(weight[i], decay[i], 0, 3, 0, m_position);
        cg_zzz->add_primitive(weight[i], decay[i], 0, 0, 3, m_position);
        cg_xxy->add_primitive(weight[i], decay[i], 2, 1, 0, m_position);
        cg_xxz->add_primitive(weight[i], decay[i], 2, 0, 1, m_position);
        cg_xyy->add_primitive(weight[i], decay[i], 1, 2, 0, m_position);
        cg_yyz->add_primitive(weight[i], decay[i], 0, 2, 1, m_position);
        cg_xzz->add_primitive(weight[i], decay[i], 1, 0, 2, m_position);
        cg_yzz->add_primitive(weight[i], decay[i], 0, 1, 2, m_position);
        cg_xyz->add_primitive(weight[i], decay[i], 1, 1, 1, m_position);
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