#pragma once

#include "BasisSet/orbital.hpp"
#include "Eigen/Dense"


class Integrator {
    public:
        virtual double electron_electron_integral(const Orbital& o1, const Orbital& o2, const Orbital& o3, const Orbital& o4) const;
        virtual double overlap_integral(const Orbital& o1, const Orbital& o2) const;
        virtual double electron_nucleus_integral(const Orbital& o1, const Orbital& o2, const Eigen::Vector3d& position) const;
        virtual double laplacian_integral(const Orbital& o1, const Orbital& o2) const;
};