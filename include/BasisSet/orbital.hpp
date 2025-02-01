#pragma once

class Orbital {
    public:
        ~Orbital() = default;
        double normalization() const { return m_normalization_constant; };

    protected:
        double m_normalization_constant = 0;
};