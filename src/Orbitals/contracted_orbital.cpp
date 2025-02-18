#include "Orbitals/contracted_orbital.hpp"
#include "Orbitals/gaussian_primitive.hpp"
#include "Orbitals/slater_primitive.hpp"

template class ContractedOrbital<GaussianPrimitive>;
template class ContractedOrbital<SlaterPrimitive>;