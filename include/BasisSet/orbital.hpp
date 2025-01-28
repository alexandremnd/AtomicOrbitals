#pragma once

class Orbital {
    public:
        virtual ~Orbital() = default;
        virtual double normalization() const;

    protected:
        double m_normalization_constant = 0;
};