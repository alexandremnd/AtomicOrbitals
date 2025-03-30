#include <cmath>
#include "Eigen/Dense"

#include "Atom/atom_list.hpp"
#include "Atom/atom.hpp"
#include "Atom/molecule.hpp"

#include "HartreeFock/restricted_hartree_fock.hpp"
#include "HartreeFock/unrestricted_hartree_fock.hpp"

#include "Numerov/sol_radiale.hpp"
#include "Numerov/minimisation.hpp"

#include "Utils/clock.hpp"
#include "Utils/misc.hpp"

void print_openmp_state() {
#ifdef _OPENMP
    std::cout << "[OpenMP] Parellization enabled. Max threads: "
              << omp_get_max_threads() << std::endl;
#else
    std::cout << "[OpenMP] No parellization enabled. Consider using smaller "
                 "basis set."
              << std::endl;
#endif
}

void print_copyright() {
    std::cout << "===========================" << std::endl;
    std::cout << "AtomicOrbitals Copyright (C) 2025 Alexandre Menard, Ewan "
                 "Bataille, Jeremy Atané"
              << std::endl;
    std::cout << "This program comes with ABSOLUTELY NO WARRANTY" << std::endl;
    std::cout << "This is free software, and you are welcome to redistribute it"
                 " under certain conditions;"
              << std::endl;
    std::cout << "See the GNU General Public License v3 for details."
              << std::endl;
    std::cout << "===========================" << std::endl;
}

void compute_periodic_table_energies(int min_Z, int max_Z,
                                     std::string basis_name) {

    auto output_file = create_output_file("atom-" + basis_name + ".out");

    auto clock = Clock();
    for (int Z = min_Z; Z < max_Z + 1; Z += 2) {
        std::cout << "Computing " << get_element_long_name(Element(Z))
                  << " ground state." << std::endl;

        auto atom = Atom(Element(Z), basis_name);

        auto rhf = RestrictedHartreeFock(atom, Z);
        rhf.set_smoothing_factor(0.7);
        rhf.run(1e-4, 1000, 2);

        *output_file << Z << "," << rhf.get_final_energy() << std::endl;
    }
    clock.time_s("Periodic table computation time: ");
}

void compute_atoms() {
    compute_periodic_table_energies(2, 16, "sto6g");
    compute_periodic_table_energies(2, 16, "6-31g");

    execute_python_script("plot_atom.py");
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
            *output_stream << step * i + min << "," << rhf.get_final_energy()
                           << std::endl;
        }
    }

    auto H2 = Molecule();
    auto Ha = std::make_shared<Atom>(Element::H, "6-31g");
    auto Hb = std::make_shared<Atom>(Element::H, "6-31g");

    H2.add_atom(Ha);
    H2.add_atom(Hb);

    std::string output = "h2_rhf.out";
    auto output_stream = create_output_file(output);

    for (int i = 0; i < n; i++) {
        Hb->set_position({0.0, 0.0, i * step + min});
        auto rhf = RestrictedHartreeFock(H2, 2);
        rhf.set_smoothing_factor(0.7);
        rhf.set_silent(true);
        rhf.run(1e-6, 1000, 1);
        *output_stream << step * i + min << "," << rhf.get_final_energy()
                       << std::endl;
    }

    execute_python_script("plot_optimization.py");
}

void compute_h2_eq() {
    auto H2 = Molecule();
    auto Ha = std::make_shared<Atom>(Element::H, "6-31g");
    auto Hb = std::make_shared<Atom>(Element::H, "6-31g");

    H2.add_atom(Ha);
    H2.add_atom(Hb);

    Hb->set_position({0.0, 0.0, 1.5});

    auto rhf = RestrictedHartreeFock(H2, 2);
    rhf.set_smoothing_factor(0.7);
    rhf.run(1e-6, 1000, 1);

    auto system_file = create_output_file("h2/system.out");
    auto density_file = create_output_file("h2/density.out");
    auto energy_file = create_output_file("h2/energy.out");
    auto coefficients_file = create_output_file("h2/coefficients.out");

    *system_file << H2 << std::endl;
    *density_file << rhf.density_matrix() << std::endl;
    *energy_file << rhf.orbital_energies() << std::endl;
    *coefficients_file << rhf.coefficient_matrix() << std::endl;

    std::cout << "Data saved in h2 directory." << std::endl;
    std::cout << "===========================" << std::endl;
    execute_python_script("plot_molecule.py");
}

/**
 * @brief Compute ethylen orbital energies and density matrix at equilibrium
 * geometry. 121.3° bond angle, 109 pm (2.06 a0) C-H bond length, 134 pm (2.53
 * a0) C=C bond length (alkene).
 *
 */
void compute_ethylen_eq() {
    auto ethylen = Molecule();
    auto C1 = std::make_shared<Atom>(Element::C, "6-311++(2d,2p)");
    auto C2 = std::make_shared<Atom>(Element::C, "6-311++(2d,2p)");
    auto H1 = std::make_shared<Atom>(Element::H, "6-311++g");
    auto H2 = std::make_shared<Atom>(Element::H, "6-311++g");
    auto H3 = std::make_shared<Atom>(Element::H, "6-311++g");
    auto H4 = std::make_shared<Atom>(Element::H, "6-311++g");

    double cos_60 = cos(0.5 * 121.3 * M_PI / 180);
    double sin_60 = sin(0.5 * 121.3 * M_PI / 180);

    C1->set_position({0, 0, 0});
    C2->set_position({0, 0, 2.53});

    // Hydrogen linked to C2
    H3->set_position({0, 2.06 * sin_60, 2.53 + 2.06 * cos_60});
    H4->set_position({0, -2.06 * sin_60, 2.53 + 2.06 * cos_60});

    // Hydrogen linked to C1
    H1->set_position({0, 2.06 * sin_60, -2.06 * cos_60});
    H2->set_position({0, -2.06 * sin_60, -2.06 * cos_60});

    ethylen.add_atom(C1);
    ethylen.add_atom(C2);
    ethylen.add_atom(H1);
    ethylen.add_atom(H2);
    ethylen.add_atom(H3);
    ethylen.add_atom(H4);

    auto rhf = RestrictedHartreeFock(ethylen, 16);
    rhf.set_silent(false);
    rhf.set_smoothing_factor(0.7);
    rhf.run(1e-6, 1000, 1);

    // Print results
    auto system_file = create_output_file("ethylen/system.out");
    auto density_file = create_output_file("ethylen/density.out");
    auto energy_file = create_output_file("ethylen/energy.out");
    auto coefficients_file = create_output_file("ethylen/coefficients.out");

    *system_file << ethylen << std::endl;
    *density_file << rhf.density_matrix() << std::endl;
    *energy_file << rhf.orbital_energies() << std::endl;
    *coefficients_file << rhf.coefficient_matrix() << std::endl;

    std::cout << "Data saved in ethylen directory." << std::endl;
    std::cout << "===========================" << std::endl;

    execute_python_script("plot_molecule.py");
}

void numerov() {
    minimisation(2, 1);
    radial_solution(2, 1, 2);

    execute_python_script("affichage.py");
}

int main() {
    print_copyright();
    print_openmp_state();

    int choice = 0;
    bool running = true;

    while (running) {
        // Display menu
        std::cout << "\n===== AtomicOrbitals Menu =====\n";
        std::cout << "1. Compute H2 at equilibrium geometry\n";
        std::cout << "2. Compute ethylene at equilibrium geometry\n";
        std::cout << "3. Run Numerov's method for hydrogen-like atoms\n";
        std::cout << "4. Compute periodic table energies\n";
        std::cout << "5. Optimize H2 bond length\n";
        std::cout << "0. Exit\n";
        std::cout << "Enter your choice: ";

        // Get user input with validation
        if (!(std::cin >> choice)) {
            std::cin.clear(); // Clear error flags
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(),
                            '\n'); // Discard invalid input
            std::cout << "Invalid input. Please enter a number.\n";
            continue;
        }

        // Process user choice
        switch (choice) {
        case 1:
            std::cout << "\nComputing H2 at equilibrium geometry...\n";
            compute_h2_eq();
            return 0;

        case 2:
            std::cout << "\nComputing ethylene at equilibrium geometry...\n";
            compute_ethylen_eq();
            return 0;

        case 3:
            std::cout
                << "\nRunning Numerov's method for hydrogen-like atoms...\n";
            numerov();
            return 0;

        case 4: {
            std::cout << "\nComputing periodic table energies...\n";

            compute_atoms();
            return 0;
        }

        case 5: {
            double min = 0.5, max = 10.0, step = 0.05;
            std::vector<std::string> basis_sets = {"6-31g", "sto6g", "3-21g"};

            optimize_h2(min, max, step, basis_sets);
            return 0;
        }

        case 0:
            std::cout << "Exiting program. Goodbye!\n";
            running = false;
            break;

        default:
            std::cout << "Invalid choice. Please try again.\n";
            break;
        }
    }

    return 0;
};