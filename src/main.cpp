#include "Atom/atom_list.hpp"
#include "Eigen/Dense"
#include "HartreeFock/restricted_hartree_fock.hpp"
#include "HartreeFock/unrestricted_hartree_fock.hpp"
#include "Atom/atom.hpp"
#include "Atom/molecule.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include "Utils/clock.hpp"

namespace fs = std::filesystem;

std::ostream create_output_file(std::string output, int precision = 10) {
    fs::path path = fs::current_path() / "data" / output;
    fs::create_directories(path.parent_path());

    std::ofstream output_file(path);
    if (!output_file.is_open()) {
        std::cerr << "Ouput file creation: Could not open the output file "
                  << path << ".\n";
        exit(1);
    }
    output_file << std::setprecision(precision);
}

void optimize_h2(double min, double max, double step,
                 std::vector<std::string> basis_sets) {

    int n = (max - min) / step;

    for (auto &basis_name : basis_sets) {
        std::cout << "Optimizing H2 with basis set " << basis_name << std::endl;

        auto H2 = Molecule();
        auto Ha = std::make_shared<Atom>(Element::H, basis_name);
        auto Hb = std::make_shared<Atom>(Element::H, basis_name);

        H2.add_atom(Ha);
        H2.add_atom(Hb);

        std::string output = "h2_" + basis_name + ".out";
        auto output_stream = create_output_file(output);

        for (int i = 0; i < n; i++) {
            Hb->set_position({0.0, 0.0, i * step + min});
            auto rhf = UnrestrictedHartreeFock(H2, 2, 1);
            rhf.set_smoothing_factor(0.7);
            rhf.set_silent(true);
            rhf.run(1e-6, 1000, 1);
            output_stream << step * i + min << "," << rhf.get_final_energy()
                          << std::endl;
        }
    }
}

void validation_table() {
    auto H2 = Molecule();
    auto Ha = std::make_shared<Atom>(Element::H, "3-21g");
    auto Hb = std::make_shared<Atom>(Element::H, "3-21g");
}

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

        auto atom = Atom(Element(Z), basis_name);

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

    // compute_periodic_table_energies(2, 8, "sto6g");*

    return 0;
};