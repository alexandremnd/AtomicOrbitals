#pragma once

class Orbital {
    public:
        ~Orbital() = default;
        double normalization() const;

    protected:
        double m_normalization_constant = 0;
};