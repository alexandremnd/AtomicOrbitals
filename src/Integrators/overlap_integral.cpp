#include "Integrators/overlap_integral.hpp"
#include <stdexcept>
#include <cassert>

double overlap_integral(const GaussianPrimitive& orbital1, const GaussianPrimitive& orbital2) {
    throw std::logic_error("Not implemented");
}

double overlap_integral(const ContractedGaussian& orbital1, const ContractedGaussian& orbital2) {
    throw std::logic_error("Not implemented");
}

double overlap_integral(const SlaterPrimitive& orbital1, const SlaterPrimitive& orbital2, const int n_offset) {
    assert(orbital1.n() + orbital2.n() + n_offset >= 0);

    if (orbital1.l() != orbital2.l() || orbital1.m() != orbital2.m()) {
        return 0.;
    }

    return boost::math::tgamma(orbital1.n() + orbital2.n() + n_offset + 1) * std::pow(orbital1.alpha() + orbital2.alpha(), -orbital1.n() - orbital2.n() - n_offset - 1);
}