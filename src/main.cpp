#include "Atom/atom.hpp"
#include "BasisSet/contracted_orbital.interface.hpp"
#include "BasisSet/slater_primitive.hpp"
#include "Eigen/Dense"
#include "HartreeFock/hamiltonian.hpp"
#include "HartreeFock/restricted_hartree_fock.hpp"

int main() {
    // auto He = Atom<SlaterPrimitive>(2, Eigen::Vector3d(0, 0, 0));
    // He.add_orbital(1, 0, 0, 2);
    // He.add_orbital(1, 0, 0, 1.6);
    // He.add_orbital(1, 0, 0, 0.3);

    // auto He2 = Atom<ContractedSlater>(2, Eigen::Vector3d(0, 0, 0));
    // He2.add_contracted_slater({1.3479, -0.001613, -0.100506, -0.270779},
    //                           {1, 3, 2, 2}, {0, 0, 0, 0}, {0, 0, 0, 0},
    //                           {1.4595, 5.3244, 2.6298, 1.7504});

    // auto Be = Atom<ContractedSlater>(4, Eigen::Vector3d(0, 0, 0));
    // Be.add_contracted_slater(
    //     {0.285107, 0.474813, -0.001620, 0.052852, 0.243499, 0.000106,
    //      -0.000032},
    //     {1, 1, 3, 3, 2, 2, 2}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0},
    //     {5.7531, 3.7156, 9.9670, 3.7128, 4.4661, 1.2919, 0.8555});
    // Be.add_contracted_slater(
    //     {-0.016378, -0.155066, 0.000426, -0.059234, -0.031925, 0.387968,
    //      0.685674},
    //     {1, 1, 3, 3, 2, 2, 2}, {0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0},
    //     {5.7531, 3.7156, 9.9670, 3.7128, 4.4661, 1.2919, 0.8555});

    // Hamiltonian hamiltonian(Be);
    // hamiltonian.print();

    // auto Be2 = Atom<SlaterPrimitive>(4, Eigen::Vector3d(0, 0, 0));
    // Be2.add_orbital(1, 0, 0, 3.7156);
    // Be2.add_orbital(2, 0, 0, 2.74);

    // RestrictedHartreeFock hf(He2, 2);

    // hf.set_smoothing_factor(0.3);
    // hf.set_diis_size(0);
    // hf.run(1e-5, 1000, 1, false);

    auto gp1 = GaussianPrimitive(0, 0, 0, 1.00, Eigen::Vector3d{1, 0, 0});
    auto gp2 = GaussianPrimitive(0, 0, 0, 1.00, Eigen::Vector3d{1, 0, 0});

    std::cout << overlap_integral(gp1, gp1) << std::endl;

    return 0;
};