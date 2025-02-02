#include "BasisSet/contracted_orbital.hpp"
#include "BasisSet/gaussian_primitive.hpp"
#include "BasisSet/slater_primitive.hpp"

template class ContractedOrbital<GaussianPrimitive>;
template class ContractedOrbital<SlaterPrimitive>;