#pragma once

#include <concepts>
#include "BasisSet/slater_primitive.hpp"
#include "BasisSet/gaussian_primitive.hpp"
#include "BasisSet/contracted_orbital.hpp"
#include "BasisSet/orbital.hpp"

// Declare the template class for all orbital type
#define DECLARE_TEMPLATE(TemplateClass) \
    template class TemplateClass<SlaterPrimitive>; \
    template class TemplateClass<ContractedSlater>; \
    // template class TemplateClass<ContractedGaussian>;

#define DECLARE_EXTERN_TEMPLATE(TemplateClass) \
    extern template class TemplateClass<SlaterPrimitive>; \
    extern template class TemplateClass<ContractedSlater>; \
    // extern template class TemplateClass<ContractedGaussian>;

// Declare the template class for the primitive types.
#define DECLARE_EXTERN_PRIMITIVE(TemplateClass) \
    extern template class TemplateClass<SlaterPrimitive>; \
    extern template class TemplateClass<GaussianPrimitive>;

#define DECLARE_TEMPLATE_PRIMITIVE(TemplateClass) \
    template class TemplateClass<SlaterPrimitive>; \
    template class TemplateClass<GaussianPrimitive>;

template <typename T>
concept DerivedFromOrbital = std::derived_from<T, Orbital>;