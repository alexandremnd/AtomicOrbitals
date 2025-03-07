#include "Atom/atom_list.hpp"
#include "Atom/atom.hpp"
#include "Atom/molecule.hpp"
#include "Eigen/Dense"
#include "HartreeFock/restricted_hartree_fock.hpp"
#include "Orbitals/contracted_orbital.interface.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include "Utils/clock.hpp"

namespace fs = std::filesystem;

void compute_periodic_table_energies(
    int min_Z, int max_Z, std::string basis_name,
    std::string output = "atomic_energies.out") {

    fs::path path = fs::current_path() / "data" / output;
    fs::create_directories(path.parent_path());

    std::ofstream output_file(path);
    if (!output_file.is_open()) {
        std::cerr << "Computation: Could not open the output file " << path
                  << ".\n";
        exit(1);
    }
    output_file << std::setprecision(10);

    auto clock = Clock();
    for (int Z = min_Z; Z < max_Z + 1; Z += 2) {
        std::cout << "Computing " << get_element_long_name(Element(Z))
                  << " ground state." << std::endl;

        auto atom = Atom<CGTO>(Element(Z), basis_name);
        std::cout << "Basis size: " << atom.get_orbitals().size() << std::endl;
        for (auto &co : atom.get_orbitals()) {
            std::cout << co << std::endl;
        }

        auto rhf = RestrictedHartreeFock(atom, Z);
        rhf.set_smoothing_factor(0.7);
        rhf.run(1e-4, 1000, 2);

        output_file << Z << "," << rhf.get_final_energy() << std::endl;
    }
    clock.time_s("Periodic table computation time: ");
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

    auto H2 = Molecule<CGTO>();
    auto Ha = std::make_shared<Atom<CGTO>>(Element::H, "6-31g");
    auto Hb = std::make_shared<Atom<CGTO>>(Element::H, "6-31g");

    H2.add_atom(Ha);
    H2.add_atom(Hb);

    std::string output = "molecular_optimization.out";

    fs::path path = fs::current_path() / "data" / output;
    fs::create_directories(path.parent_path());

    std::ofstream output_file(path);
    if (!output_file.is_open()) {
        std::cerr << "Computation: Could not open the output file " << path
                  << ".\n";
        exit(1);
    }
    output_file << std::setprecision(10);

    double max = 7;
    double min = 0.01;
    double step = 0.01;
    int n = (max - min) / step;

    for (int i = 0; i < n; i++) {
        Hb->set_position({0.0, 0.0, i * step + min});
        auto rhf = RestrictedHartreeFock(H2, 2);
        rhf.set_smoothing_factor(0.7);
        rhf.run(1e-6, 1000, 5);
        output_file << step * i + min << "," << rhf.get_final_energy()
                    << std::endl;
    }

    return 0;
};