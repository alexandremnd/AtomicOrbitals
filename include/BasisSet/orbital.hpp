#pragma once

class Orbital {
    public:
        Orbital();
        ~Orbital();

        virtual double normalization() const;

    protected:
        double m_normalization_constant = 0;
};