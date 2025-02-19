#pragma once

#include "Atom/atom.interface.hpp"
#include "Atom/atom_list.hpp"
#include "Orbitals/contracted_orbital.interface.hpp"
#include "Orbitals/gaussian_primitive.hpp"
#include "Orbitals/slater_primitive.hpp"

/**
 * @brief Loads a basis file and parses the basis set into the passed atom.
 * @note Basis file is expected to be in the data/BasisSet/basis_name.basis
 * directory from the current working directory. Working directory is the
 * directory from which the executable is run.
 *
 * @param elt Element to parse basis for
 * @param basis_name Name of the basis set to load
 * @param atom Atom to load the basis set into
 */
void parse_basis(Element elt, std::string basis_name,
                 Atom<ContractedGaussian> &atom);

/**
 * @brief Loads a basis file and parses the basis set into the passed atom.
 * @note Basis file is expected to be in the data/BasisSet/basis_name directory
 * from the current working directory. Working directory is the directory from
 * which the executable is run.
 *
 * @param elt Element to parse basis for
 * @param basis_name Name of the basis set to load
 * @param atom Atom to load the basis set into
 */
void parse_basis(Element elt, std::string basis_name,
                 Atom<ContractedSlater> &atom);