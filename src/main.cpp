#include "Atom/atom_list.hpp"
#include "Atom/atom.hpp"
#include "Eigen/Dense"
#include "HartreeFock/hamiltonian.hpp"
#include "HartreeFock/restricted_hartree_fock.hpp"
#include "Orbitals/contracted_orbital.interface.hpp"
#include "Orbitals/gaussian_primitive.hpp"
#include "Utils/clock.hpp"
#include <omp.h>

int main() {
    // auto He = Atom<SlaterPrimitive>(2, Eigen::Vector3d(0, 0, 0));
    // He.add_orbital(1, 0, 0, 2);
    // He.add_orbital(1, 0, 0, 1.6);
    // He.add_orbital(1, 0, 0, 0.3);

    // auto He2 = Atom<ContractedSlater>(2, Eigen::Vector3d(0, 0, 0));
    // He2.add_contracted_slater({1.3479, -0.001613, -0.100506, -0.270779},
    //                           {1, 3, 2, 2}, {0, 0, 0, 0}, {0, 0, 0, 0},
    //                           {1.4595, 5.3244, 2.6298, 1.7504});

    auto atom = Atom<CGTO>(Element(24), "sto6g");

    auto rhf = RestrictedHartreeFock(atom, 24);
    rhf.set_smoothing_factor(0.9);
    rhf.run(1e-6, 1000, 2, false);
    // clock.time_ms("Execution time");

    return 0;
};