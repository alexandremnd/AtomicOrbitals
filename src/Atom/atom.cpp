#include "include/Atom/atom.hpp"

float Atom::getEnergyLevel(int n)
{
    return - (m_Z * m_Z) / (2. * n * n);
}