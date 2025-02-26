#include "Atom/atom_list.hpp"
#include "Orbitals/contracted_orbital.interface.hpp"
#include "Utils/string_fun.hpp"
#include "Atom/atom.hpp"

#include <fstream>
#include <ostream>
#include <iostream>
#include <string>
#include <filesystem>

namespace fs = std::filesystem;

template <>
void parse_basis(Element elt, std::string basis_name,
                 Atom<ContractedGaussian> &atom) {
    fs::path basis_path = fs::current_path() / "data" / "BasisSet" /
                          (to_lowercase(basis_name) + ".basis");

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

    bool found = false;
    std::string atom_to_find = get_element_short_name(elt);

    std::string orbital_type;
    int orbital_size = 0;
    std::vector<double> weight, weight_sp, decay;

    std::string line;
    while (std::getline(basis_file, line)) {
        std::istringstream iss(line);

        // Skip non relevant lines until we find the atom we are looking for.
        if (!found) {
            std::string possible_atom_name;
            iss >> possible_atom_name;

            if (possible_atom_name == atom_to_find) {
                found = true;
            }

            continue;
        }

        // We arrived at the end of the basis set for the desired atom
        if (line == "****" && found) {
            break;
        }

        // We found the atom we are looking for, but we are not at the end of
        // the basis set. Parse the basis set.

        std::string first_word, second_word, third_word;
        iss >> first_word;
        iss >> second_word;
        iss >> third_word;

        if (first_word == "S" || first_word == "P" || first_word == "D" ||
            first_word == "F" || first_word == "SP") {

            // Add the currently built orbital to the atom
            if (orbital_type == "S") {
                atom.add_gaussian_orbital_stype(weight, decay);
            } else if (orbital_type == "P") {
                atom.add_gaussian_orbital_ptype(weight, decay);
            } else if (orbital_type == "D") {
                atom.add_gaussian_orbital_dtype(weight, decay);
            } else if (orbital_type == "F") {
                atom.add_gaussian_orbital_ftype(weight, decay);
            } else if (orbital_type == "SP") {
                atom.add_gaussian_orbital_stype(weight, decay);
                atom.add_gaussian_orbital_ptype(weight_sp, decay);
            }

            // Reset the type, weight, and decay vectors for next orbital
            orbital_type = first_word;
            iss >> orbital_size;

            weight.clear();
            weight_sp.clear();
            decay.clear();

            weight.reserve(N);
            weight_sp.reserve(N);
            decay.reserve(N);

            continue;
        }

        replace_letter(first_word, 'D', 'E');
        replace_letter(second_word, 'D', 'E');

        double alpha = std::stod(first_word);
        double c = std::stod(second_word);

        weight.push_back(c);
        decay.push_back(alpha);

        if (orbital_type == "SP") {
            replace_letter(third_word, 'D', 'E');
            double c_sp = std::stod(third_word);
            weight_sp.push_back(c_sp);
        }
    }

    // Add the last orbital to the atom
    if (orbital_type == "S") {
        atom.add_gaussian_orbital_stype(weight, decay);
    } else if (orbital_type == "P") {
        atom.add_gaussian_orbital_ptype(weight, decay);
    } else if (orbital_type == "D") {
        atom.add_gaussian_orbital_dtype(weight, decay);
    } else if (orbital_type == "F") {
        atom.add_gaussian_orbital_ftype(weight, decay);
    } else if (orbital_type == "SP") {
        atom.add_gaussian_orbital_stype(weight, decay);
        atom.add_gaussian_orbital_ptype(weight, decay);
    }

    return;
}

template <>
void parse_basis(Element elt, std::string basis_name,
                 Atom<ContractedSlater> &atom) {

    return;
}