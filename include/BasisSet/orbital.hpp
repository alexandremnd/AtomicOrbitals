#pragma once

class Orbital {
    public:
        virtual ~Orbital() = default;
        virtual double normalization() const;

        // virtual double overlap(const Orbital& other) const;
        // virtual double laplacian(const Orbital& other) const;
        // virtual double electron_nucleus(const Orbital& other, const Eigen::Vector3d& position) const;
        // virtual double electron_electron(const Orbital& o2, const Orbital& o3, const Orbital& o4) const;

    protected:
        double m_normalization_constant = 0;
};