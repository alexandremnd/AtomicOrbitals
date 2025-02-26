#include "Atom/atom_list.hpp"
#include "Atom/atom.hpp"
#include "Eigen/Dense"
#include "HartreeFock/hamiltonian.hpp"
#include "HartreeFock/restricted_hartree_fock.hpp"
#include "Orbitals/contracted_orbital.interface.hpp"
#include "Orbitals/gaussian_primitive.hpp"
#include "Utils/clock.hpp"
#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;

void compute_atom_energies() {
    fs::path data_dir = fs::current_path() / "data" / "atom_energies.out";
    std::ofstream outfile(data_dir);

    for (int i = 2; i < 37; i += 2) {
        if (i == 16)
            continue;
        auto atom = Atom<CGTO>(Element(i), "3-21g");
        auto rhf = RestrictedHartreeFock(atom, i);
        rhf.set_smoothing_factor(0.9);
        rhf.run(1e-6, 1000, 2, true);

        double energy = rhf.get_final_energy();
        outfile << "Atom: " << i << " Energy: " << energy << std::endl;
    }

    outfile.close();
}

int main() {
#ifdef _OPENMP
    std::cout << "[OpenMP] Parellization enabled. Max threads: "
              << omp_get_max_threads() << std::endl;
#else
    std::cout << "[OpenMP] No parellization enabled. Consider using smaller "
                 "basis set."
              << std::endl;
#endif
    // compute_atom_energies();

    auto atom = Atom<CGTO>(Element(14), "ugbs");

    auto rhf = RestrictedHartreeFock(atom, 14);
    rhf.set_smoothing_factor(0.3);
    rhf.run(1e-6, 1000, 2, false);

    return 0;
};