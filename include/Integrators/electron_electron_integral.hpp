#pragma once

#include "BasisSet/slater_primitive.hpp"
#include "BasisSet/gaussian_primitive.hpp"
#include "BasisSet/contracted_orbital.hpp"
#include "concepts.hpp"


double electron_electron_integral(const GaussianPrimitive& orbital1,
                                const GaussianPrimitive& orbital2,
                                const GaussianPrimitive& orbital3,
                                const GaussianPrimitive& orbital4);

double electron_electron_integral(const SlaterPrimitive& orbital1,
                                const SlaterPrimitive& orbital2,
                                const SlaterPrimitive& orbital3,
                                const SlaterPrimitive& orbital4);

template <DerivedFromOrbital PrimitiveType>
double electron_electron_integral(const ContractedOrbital<PrimitiveType>& orbital1,
                                const ContractedOrbital<PrimitiveType>& orbital2,
                                const ContractedOrbital<PrimitiveType>& orbital3,
                                const ContractedOrbital<PrimitiveType>& orbital4) {
    double integral = 0.0;

    for (size_t i = 0; i < orbital1.size(); i++) {
        for (size_t j = 0; j < orbital2.size(); j++) {
            for (size_t k = 0; k < orbital3.size(); k++) {
                for (size_t l = 0; l < orbital4.size(); l++) {
                    integral += orbital1.get_coefficient(i) * orbital2.get_coefficient(j) * orbital3.get_coefficient(k) * orbital4.get_coefficient(l) *
                                electron_electron_integral(orbital1.get_primitive(i), orbital2.get_primitive(j), orbital3.get_primitive(k), orbital4.get_primitive(l));
                }
            }
        }
    }

    return integral * orbital1.normalization() * orbital2.normalization() * orbital3.normalization() * orbital4.normalization();
}