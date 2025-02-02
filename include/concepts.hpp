#pragma once

#include "BasisSet/gaussian_primitive.hpp"
#include "BasisSet/orbital.hpp"
#include "BasisSet/slater_primitive.hpp"
#include <concepts>

// Declare the template class for all orbital type
#define DECLARE_TEMPLATE(TemplateClass)                                        \
    template class TemplateClass<SlaterPrimitive>;                             \
    template class TemplateClass<ContractedSlater>;                            \
    // template class TemplateClass<ContractedGaussian>;

#define DECLARE_EXTERN_TEMPLATE(TemplateClass)                                 \
    extern template class TemplateClass<SlaterPrimitive>;                      \
    extern template class TemplateClass<ContractedSlater>;                     \
    // extern template class TemplateClass<ContractedGaussian>;

template <typename T>
concept DerivedFromOrbital = std::derived_from<T, Orbital>;