#pragma once

#include "Atom/atom.hpp"
#include "Atom/atom_list.hpp"

/**
 * @brief Loads a basis file and parses the basis set into the passed atom.
 * @note Basis file is expected to be in the data/BasisSet/basis_name.basis
 * directory from the current working directory. Working directory is the
 * directory from which the executable is run.
 *
 * @note The basis file is expected to be in the Gaussian basis set format.
 *
 * @param elt Element to parse basis for
 * @param basis_name Name of the basis set to load
 * @param atom Atom to load the basis set into
 */
void parse_basis(Element elt, std::string basis_name, Atom &atom);
