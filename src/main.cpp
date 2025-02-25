#include "Atom/atom_list.hpp"
#include "Atom/atom.hpp"
#include "Eigen/Dense"
#include "HartreeFock/hamiltonian.hpp"
#include "HartreeFock/restricted_hartree_fock.hpp"
#include "Orbitals/contracted_orbital.interface.hpp"
#include "Orbitals/gaussian_primitive.hpp"
#include "Utils/clock.hpp"

int main() {
#ifdef _OPENMP
    std::cout << "[OpenMP] Parellization enabled. Max threads: "
              << omp_get_max_threads() << std::endl;
#else
    std::cout << "[OpenMP] No parellization enabled. Consider using smaller "
                 "basis set."
              << std::endl;
#endif

    auto atom = Atom<CGTO>(Element(10), "3-21g");

    for (auto &orb : atom.get_orbitals()) {
        std::cout << orb << std::endl;
    }

    auto rhf = RestrictedHartreeFock(atom, 10);
    rhf.set_smoothing_factor(0.9);
    rhf.run(1e-6, 1000, 2, false);
    // clock.time_ms("Execution time");

    return 0;
};