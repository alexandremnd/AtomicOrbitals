#pragma once

#include "Integrators/overlap_primitive.hpp"
#include "BasisSet/contracted_orbital.hpp"
#include "concepts.hpp"

template <DerivedFromOrbital PrimitiveType>
double overlap_integral(const ContractedOrbital<PrimitiveType>& orbital1, const ContractedOrbital<PrimitiveType>& orbital2) {
    double integral = 0.0;

    for (size_t i = 0; i < orbital1.size(); i++) {
        for (size_t j = 0; j < orbital2.size(); j++) {
            integral += orbital1.get_coefficient(i) * orbital2.get_coefficient(j) *
                        overlap_integral(orbital1.get_primitive(i), orbital2.get_primitive(j));
        }
    }

    return integral * orbital1.normalization() * orbital2.normalization();
}