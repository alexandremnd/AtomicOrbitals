import sympy as sp
import numpy as np


def symbolic_overlap(n1, l1, m1, a1, n2, l2, m2, a2):
    """
    Calculate the overlap between two hydrogenic orbitals.

    Parameters
    ----------
    n1 : int
        Principal quantum number of the first orbital.
    l1 : int
        Azimuthal quantum number of the first orbital.
    m1 : int
        Magnetic quantum number of the first orbital.
    a1 : float
        Exponential decay of the first orbital.
    n2 : int
        Principal quantum number of the second orbital.
    l2 : int
        Azimuthal quantum number of the second orbital.
    m2 : int
        Magnetic quantum number of the second orbital.
    a2 : float
        Exponential decay of the second orbital.

    Returns
    -------
    float
        The overlap between the two orbitals.
    """
    r, theta, phi = sp.symbols('r theta phi')
    Y1 = sp.Ynm(l1, m1, theta, phi).expand(func=True)
    Y2 = sp.Ynm(l2, m2, theta, phi).expand(func=True)

    return sp.integrate(r**(n1 + n2) * sp.exp(-(a1 + a2) * r) * sp.conjugate(Y1) * Y2 * sp.sin(theta), (r, 0, sp.oo), (theta, 0, np.pi), (phi, 0, 2 * np.pi))

if __name__ == "__main__":
    print(symbolic_overlap(1, 0, 0, 0.5, 5, 0, 0, 1.3).n())
    print(symbolic_overlap(3, 0, 0, 0.5, 8, 0, 0, 1.3).n())