#include "parse_basis.hpp"
#include "Atom/atom_list.hpp"
#include "Orbitals/contracted_orbital.interface.hpp"
#include "Utils/to_lowercase.hpp"
#include "Atom/atom.hpp"

#include <fstream>
#include <ostream>
#include <iostream>
#include <string>
#include <filesystem>

namespace fs = std::filesystem;

void parse_basis(Element elt, std::string basis_name,
                 Atom<ContractedGaussian> &atom) {
    fs::path basis_path = fs::current_path() / "data" / "BasisSet" /
                          to_lowercase(basis_name) /
                          (to_lowercase(get_element_name(elt)) + ".basis");

    if (!fs::exists(basis_path)) {
        std::cerr << "parse_basis: The basis file " << basis_path
                  << " does not exist.\n";
        exit(1);
    }

    std::ifstream basis_file(basis_path);
    if (!basis_file.is_open()) {
        std::cerr << "parse_basis: Could not open the basis file " << basis_path
                  << ".\n";
        exit(1);
    }

    std::string line;
    while (std::getline(basis_file, line)) {
        std::istringstream iss(line);
        std::string type;
        int N;

        iss >> type;
        iss >> N;

        std::vector<double> weight(N), decay(N);

        for (int i = 0; i < N; i++) {
            iss >> weight[i];
        }

        for (int i = 0; i < N; i++) {
            iss >> decay[i];
        }

        if (type == "S") {
            atom.add_gaussian_orbital_stype(weight, decay);
        } else if (type == "P") {
            atom.add_gaussian_orbital_ptype(weight, decay);
        } else if (type == "D") {
            atom.add_gaussian_orbital_dtype(weight, decay);
        } else if (type == "F") {
            atom.add_gaussian_orbital_ftype(weight, decay);
        }
    }

    return;
}

void parse_basis(Element elt, std::string basis_name,
                 Atom<ContractedSlater> &atom) {
    // fs::path basis_path = fs::current_path() / "data" / "BasisSet" /
    //                       to_lowercase(basis_name) /
    //                       (to_lowercase(get_element_name(elt)) + ".basis");

    // if (!fs::exists(basis_path)) {
    //     std::cerr << "parse_basis: The basis file " << basis_path
    //               << " does not exist.\n";
    //     exit(1);
    // }

    // std::ifstream basis_file(basis_path);
    // if (!basis_file.is_open()) {
    //     std::cerr << "parse_basis: Could not open the basis file " <<
    //     basis_path
    //               << ".\n";
    //     exit(1);
    // }

    // std::string line;
    // while (std::getline(basis_file, line)) {
    //     std::istringstream iss(line);
    //     std::string type;
    //     int n, N; // N is the number of orbtitals, n is the principal quantum
    //               // number

    //     iss >> n;
    //     iss >> type;
    //     iss >> N;

    //     std::vector<double> weight(N), decay(N);

    //     for (int i = 0; i < N; i++) {
    //         iss >> weight[i];
    //     }

    //     for (int i = 0; i < N; i++) {
    //         iss >> decay[i];
    //     }

    //     atom.add_slater_orbital(weight, decay);
    //     atom.add_slater_orbital(weight, decay);
    //     atom.add_slater_orbital(weight, decay);
    //     atom.add_slater_orbital(weight, decay);
    // }

    return;
}