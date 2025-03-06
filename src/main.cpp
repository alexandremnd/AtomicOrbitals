#include "Atom/atom_list.hpp"
#include "Atom/atom.hpp"
#include "Eigen/Dense"
#include "HartreeFock/restricted_hartree_fock.hpp"
#include "Orbitals/contracted_orbital.interface.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>

namespace fs = std::filesystem;

void compute_periodic_table_energies(
    int min_Z, int max_Z, std::string basis_name,
    std::string output = "atomic_energies.out") {

    fs::path path = fs::current_path() / "data" / output;

    if (!fs::exists(path)) {
        std::cerr << "parse_basis: The basis file " << path
                  << " does not exist.\n";
        exit(1);
    }

    std::fstream basis_file(path);
    if (!basis_file.is_open()) {
        std::cerr << "parse_basis: Could not open the basis file " << path
                  << ".\n";
        exit(1);
    }
    basis_file << std::setprecision(10);

    for (int Z = min_Z; Z < max_Z + 1; Z += 2) {
        std::cout << "Computing " << get_element_long_name(Element(Z))
                  << " ground state." << std::endl;

        auto atom = Atom<CGTO>(Element(Z), basis_name);

        auto rhf = RestrictedHartreeFock(atom, Z);
        rhf.set_smoothing_factor(0.7);
        rhf.run(1e-4, 1000, 2);

        basis_file << Z << "," << rhf.get_final_energy() << std::endl;
    }
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

    compute_periodic_table_energies(2, 36, "ugbs");

    return 0;
};