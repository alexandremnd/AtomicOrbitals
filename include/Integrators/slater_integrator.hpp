#pragma once

#include "Integrators/integrator.hpp"

class SlaterIntegrator : public Integrator {
    public:
        double electron_electron_integral(const Orbital& o1, const Orbital& o2, const Orbital& o3, const Orbital& o4) const override;
        double overlap_integral(const Orbital& o1, const Orbital& o2) const override;
        double electron_nucleus_integral(const Orbital& o1, const Orbital& o2, const Eigen::Vector3d& position) const override;
        double laplacian_integral(const Orbital& o1, const Orbital& o2) const override;
};