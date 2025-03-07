#pragma once

#include "Atom/atom.interface.hpp"
#include "Atom/atom_list.hpp"
#include "Orbitals/contracted_orbital.interface.hpp"
#include "concepts.hpp"
#include <stdexcept>

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
template <DerivedFromOrbital OrbitalType>
void parse_basis(Element elt, std::string basis_name, Atom<OrbitalType> &atom) {
    throw std::logic_error("parse_basis: Basis parsing not implemented for the "
                           "given orbital type.");
}

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
template <>
void parse_basis(Element elt, std::string basis_name,
                 Atom<ContractedGaussian> &atom);
/**
 * @brief Loads a basis file and parses the basis set into the passed atom.
 * @note Basis file is expected to be in the data/BasisSet/basis_name.basis
 * directory from the current working directory. Working directory is the
 * directory from which the executable is run.
 *
 * @note The basis file is expected to be in the custom basis set format (see
 * data/BasisSet/sto.basis).
 *
 * @param elt Element to parse basis for
 * @param basis_name Name of the basis set to load
 * @param atom Atom to load the basis set into
 */
template <>
void parse_basis(Element elt, std::string basis_name,
                 Atom<ContractedSlater> &atom);