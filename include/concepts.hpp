#pragma once

#include <concepts>
#include "BasisSet/orbital.hpp"

#define DECLARE_TEMPLATE(TemplateClass) \
    template class TemplateClass<SlaterPrimitive>; \
    template class TemplateClass<ContractedSlater>; \
    template class TemplateClass<ContractedGaussian>;

#define DECLARE_EXTERN_TEMPLATE(TemplateClass) \
    extern template class TemplateClass<SlaterPrimitive>; \
    extern template class TemplateClass<ContractedSlater>; \
    extern template class TemplateClass<ContractedGaussian>;

template <typename T>
concept DerivedFromOrbital = std::derived_from<T, Orbital>;