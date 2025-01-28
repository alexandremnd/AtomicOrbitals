#pragma once

#include "Atom/atom.interface.hpp"
#include "BasisSet/gaussian_contracted.hpp"
#include "BasisSet/slater_primitive.hpp"
#include "BasisSet/slater_contracted.hpp"
#include "concepts.hpp"

DECLARE_EXTERN_TEMPLATE(Atom)

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::add_orbital(const OrbitalType& orbital) {
    m_orbitals.push_back(std::make_unique<OrbitalType>(orbital));
}

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::add_slater_orbital(const int n, const int l, const int m, const double alpha) requires std::is_same_v<OrbitalType, SlaterPrimitive>{
    m_orbitals.push_back(std::make_unique<SlaterPrimitive>(n, l, m, alpha));
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

    std::vector<SlaterPrimitive> primitives;
    for (size_t i = 0; i < weight.size(); i++) {
        primitives.push_back(SlaterPrimitive(n[i], l[i], m[i], decay[i]));
    }

    m_orbitals.push_back(std::make_unique<ContractedSlater>(weight, primitives));
}

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::add_gaussian_orbital_stype(const std::vector<double>& weight, const std::vector<double>& decay) requires std::is_same_v<OrbitalType, ContractedGaussian> {
    if (weight.size() != decay.size()) {
        throw std::invalid_argument("Atom: The weight and decay vectors must have the same size.");
    }

    ContractedGaussian cg{};
    for (size_t i = 0; i < weight.size(); i++) {
        cg.add_primitive(weight[i], decay[i], 0, 0, 0);
    }

    m_orbitals.push_back(std::make_unique<ContractedGaussian>(cg));
}

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::add_gaussian_orbital_ptype(const std::vector<double>& weight, const std::vector<double>& decay) requires std::is_same_v<OrbitalType, ContractedGaussian> {
    if (weight.size() != decay.size()) {
        throw std::invalid_argument("Atom: The weight and decay vectors must have the same size.");
    }

    ContractedGaussian cg_x{};
    ContractedGaussian cg_y{};
    ContractedGaussian cg_z{};

    for (size_t i = 0; i < weight.size(); i++) {
        cg_x.add_primitive(weight[i], decay[i], 1, 0, 0);
        cg_y.add_primitive(weight[i], decay[i], 0, 1, 0);
        cg_z.add_primitive(weight[i], decay[i], 0, 0, 1);
    }

    m_orbitals.push_back(cg_x);
    m_orbitals.push_back(cg_y);
    m_orbitals.push_back(cg_z);
}

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::add_gaussian_orbital_dtype(const std::vector<double>& weight, const std::vector<double>& decay) requires std::is_same_v<OrbitalType, ContractedGaussian> {
    if (weight.size() != decay.size()) {
        throw std::invalid_argument("Atom: The weight and decay vectors must have the same size.");
    }

    ContractedGaussian cg_xx{};
    ContractedGaussian cg_yy{};
    ContractedGaussian cg_zz{};
    ContractedGaussian cg_xy{};
    ContractedGaussian cg_xz{};
    ContractedGaussian cg_yz{};

    for (size_t i = 0; i < weight.size(); i++) {
        cg_xx.add_primitive(weight[i], decay[i], 2, 0, 0);
        cg_yy.add_primitive(weight[i], decay[i], 0, 2, 0);
        cg_zz.add_primitive(weight[i], decay[i], 0, 0, 2);
        cg_xy.add_primitive(weight[i], decay[i], 1, 1, 0);
        cg_xz.add_primitive(weight[i], decay[i], 1, 0, 1);
        cg_yz.add_primitive(weight[i], decay[i], 0, 1, 1);
    }

    m_orbitals.push_back(cg_xx);
    m_orbitals.push_back(cg_yy);
    m_orbitals.push_back(cg_zz);
    m_orbitals.push_back(cg_xy);
    m_orbitals.push_back(cg_xz);
    m_orbitals.push_back(cg_yz);
}

template <DerivedFromOrbital OrbitalType>
void Atom<OrbitalType>::add_gaussian_orbital_ftype(const std::vector<double>& weight, const std::vector<double>& decay) requires std::is_same_v<OrbitalType, ContractedGaussian> {
    if (weight.size() != decay.size()) {
        throw std::invalid_argument("Atom: The weight and decay vectors must have the same size.");
    }

    ContractedGaussian cg_xxx{};
    ContractedGaussian cg_yyy{};
    ContractedGaussian cg_zzz{};
    ContractedGaussian cg_xxy{};
    ContractedGaussian cg_xxz{};
    ContractedGaussian cg_xyy{};
    ContractedGaussian cg_yyz{};
    ContractedGaussian cg_xzz{};
    ContractedGaussian cg_yzz{};
    ContractedGaussian cg_xyz{};

    for (size_t i = 0; i < weight.size(); i++) {
        cg_xxx.add_primitive(weight[i], decay[i], 3, 0, 0);
        cg_yyy.add_primitive(weight[i], decay[i], 0, 3, 0);
        cg_zzz.add_primitive(weight[i], decay[i], 0, 0, 3);
        cg_xxy.add_primitive(weight[i], decay[i], 2, 1, 0);
        cg_xxz.add_primitive(weight[i], decay[i], 2, 0, 1);
        cg_xyy.add_primitive(weight[i], decay[i], 1, 2, 0);
        cg_yyz.add_primitive(weight[i], decay[i], 0, 2, 1);
        cg_xzz.add_primitive(weight[i], decay[i], 1, 0, 2);
        cg_yzz.add_primitive(weight[i], decay[i], 0, 1, 2);
        cg_xyz.add_primitive(weight[i], decay[i], 1, 1, 1);
    }

    m_orbitals.push_back(cg_xxx);
    m_orbitals.push_back(cg_yyy);
    m_orbitals.push_back(cg_zzz);
    m_orbitals.push_back(cg_xxy);
    m_orbitals.push_back(cg_xxz);
    m_orbitals.push_back(cg_xyy);
    m_orbitals.push_back(cg_yyz);
    m_orbitals.push_back(cg_xzz);
    m_orbitals.push_back(cg_yzz);
    m_orbitals.push_back(cg_xyz);
}