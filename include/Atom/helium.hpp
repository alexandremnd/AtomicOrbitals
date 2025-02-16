#pragma once

#include "Atom/atom.hpp"
#include "BasisSet/contracted_orbital.hpp"

class Helium_321G : public Atom<ContractedGaussian> {
  public:
    Helium_321G() : Atom<ContractedGaussian>(2, Eigen::Vector3d(0, 0, 0)) {
        add_gaussian_orbital_stype({0.17523, 0.8934830}, {13.6267, 1.99935});
        add_gaussian_orbital_stype({1}, {0.3829930});
    }
};

class Helium_STO6G : public Atom<ContractedGaussian> {
  public:
    Helium_STO6G() : Atom<ContractedGaussian>(2, Eigen::Vector3d(0, 0, 0)) {
        add_gaussian_orbital_stype({0.00916359628, 0.04936149294, 0.1685383049,
                                    0.3705627997, 0.4164915298, 0.1303340841},
                                   {65.98456824, 12.09819836, 3.384639924,
                                    1.162715163, 0.451516322, 0.185959356});
    }
};

class Helium_STO6Gd : public Atom<ContractedGaussian> {
  public:
    Helium_STO6Gd() : Atom<ContractedGaussian>(2, Eigen::Vector3d(0, 0, 0)) {
        add_gaussian_orbital_stype({0.00916359628, 0.04936149294, 0.1685383049,
                                    0.3705627997, 0.4164915298, 0.1303340841},
                                   {65.98456824, 12.09819836, 3.384639924,
                                    1.162715163, 0.451516322, 0.185959356});
        add_gaussian_orbital_dtype({1}, {2});
    }
};

class Helium_6311pGss : public Atom<ContractedGaussian> {
  public:
    Helium_6311pGss() : Atom<ContractedGaussian>(2, Eigen::Vector3d(0, 0, 0)) {
        add_gaussian_orbital_stype({1}, {233.304073});
        add_gaussian_orbital_stype({1}, {27.47586064});
        add_gaussian_orbital_stype({1}, {5.494186309});
        add_gaussian_orbital_stype({1}, {1.390503359});
        add_gaussian_orbital_stype({1}, {0.3990010449});
        add_gaussian_orbital_stype({1}, {234.0898557});
    }
};

class Helium_6311G2df2pd : public Atom<ContractedGaussian> {
  public:
    Helium_6311G2df2pd()
        : Atom<ContractedGaussian>(2, Eigen::Vector3d(0, 0, 0)) {
        add_gaussian_orbital_stype({0.0287452, 0.2080610, 0.8376350},
                                   {98.1243000, 14.7689000, 3.3188305});
        add_gaussian_orbital_stype({1}, {0.8740470});
        add_gaussian_orbital_stype({1}, {0.2445640});

        add_gaussian_orbital_ptype({1}, {1.5});

        add_gaussian_orbital_dtype({1}, {2});
    }
};