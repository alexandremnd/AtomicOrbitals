#pragma once

#include "Atom/atom.hpp"
#include "BasisSet/contracted_orbital.hpp"

class Beryllium_STO : public Atom<SlaterPrimitive> {
  public:
    Beryllium_STO() : Atom<SlaterPrimitive>(4, Eigen::Vector3d(0, 0, 0)) {
        add_orbital(1, 0, 0, 3.7156);
        add_orbital(2, 0, 0, 2.74);
    }
};

class Beryllium_STO6G : public Atom<ContractedGaussian> {
  public:
    Beryllium_STO6G() : Atom<ContractedGaussian>(4, Eigen::Vector3d(0, 0, 0)) {
        add_gaussian_orbital_stype({0.00916359628, 0.04936149294, 0.1685383049,
                                    0.3705627997, 0.4164915298, 0.1303340841},
                                   {312.8704937, 57.36446253, 16.0485094,
                                    5.513096119, 2.140896553, 0.8817394283});
        add_gaussian_orbital_stype({-0.01325278809, -0.04699171014,
                                    -0.03378537151, 0.2502417861, 0.5951172526,
                                    0.2407061763},
                                   {13.63324744, 2.698375464, 0.8386530829,
                                    0.3226600698, 0.1401314882, 0.0642325139});

        add_gaussian_orbital_ptype({0.0037596966, 0.0376793698, 0.1738967435,
                                    0.4180364347, 0.4258595477, 0.1017082955},
                                   {13.63324744, 2.698375464, 0.8386530829,
                                    0.3226600698, 0.1401314882, 0.0642325139});
    }
};

class Beryllium_631ppGss : public Atom<ContractedGaussian> {
  public:
    Beryllium_631ppGss()
        : Atom<ContractedGaussian>(4, Eigen::Vector3d(0, 0, 0)) {

        add_gaussian_orbital_stype(
            {0.00228574, 0.0175938, 0.0863315, 0.281835, 0.640594, 0.144467},
            {1682.8, 251.715, 57.4116, 16.5171, 4.85364, 0.626863});
        add_gaussian_orbital_stype({0.108621, 0.927301, -0.00297169},
                                   {8.30938, 1.74075, 0.485816});
        add_gaussian_orbital_stype({1.0}, {0.163613});
        add_gaussian_orbital_stype({1.0}, {0.0567285});
        add_gaussian_orbital_stype({1.0}, {0.0207});

        add_gaussian_orbital_ptype({0.0361344, 0.216958, 0.841839},
                                   {8.30938, 1.74075, 0.485816});
        add_gaussian_orbital_ptype({1.0}, {0.163613});
        add_gaussian_orbital_ptype({1.0}, {0.0567285});
        add_gaussian_orbital_ptype({1.0}, {0.0207});

        add_gaussian_orbital_dtype({1.0}, {0.255});
    }
};