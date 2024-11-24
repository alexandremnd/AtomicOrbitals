#pragma once

class Atom {
    public:
        Atom(int Z) : m_Z(Z) {};

        float getEnergyLevel(int n);

        int m_Z;
};