import sympy as sp
import numpy as np

def symbolic_laplacian(n1, l1, m1, a1, n2, l2, m2, a2):
    """
    Calculate the Laplacian of the product of two hydrogenic orbitals.

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
        The Laplacian of the product of the two orbitals.
    """
    r, theta, phi = sp.symbols('r theta phi')

    Y1 = sp.conjugate(sp.Ynm(l1, m1, theta, phi).expand(func=True))
    Y2 = sp.Ynm(l2, m2, theta, phi).expand(func=True)

    R1 = r**(n1-1) * sp.exp(-a1 * r)
    R2 = r**(n2-1) * sp.exp(-a2 * r)

    psi_1 = R1 * Y1
    psi_2 = R2 * Y2
    spherical_laplacian = lambda f: 1/(r**2) * sp.diff(r**2 * sp.diff(f, r), r) + 1/(r**2 * sp.sin(theta)) * sp.diff(sp.sin(theta) * sp.diff(f, theta), theta) + 1/(r**2 * sp.sin(theta)**2) * sp.diff(f, phi, 2)

    return sp.integrate(psi_1 * spherical_laplacian(psi_2) * r**2 * sp.sin(theta),
        (r, 0, sp.oo),
        (theta, 0, np.pi),
        (phi, 0, 2 * np.pi)
    )

if __name__ == '__main__':
    print(symbolic_laplacian(1, 0, 0, 0.5, 5, 0, 0, 1.3).n())
    print(symbolic_laplacian(5, 1, 0, 0.5, 3, 1, 0, 1.3).n())
    print(symbolic_laplacian(5, 1, 1, 0.5, 3, 1, 1, 1.3).n())
    print(symbolic_laplacian(5, 2, 1, 0.5, 3, 2, 1, 1.3).n())