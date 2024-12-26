import sympy as sp
import sympy.physics.hydrogen as sh
import scipy.integrate as spi
import scipy.special as sps
import numpy as np
import warnings

def numerical_radial_integral(n1, a1, n2, a2, n3, a3, n4, a4, L):
    """Computes the radial integral of the Slater-type orbitals using numerical integration.

    Args:
        n1 (int): Principal quantum number of the first orbital.
        a1 (float): Exponential decay of the first orbital.
        n2 (int): Principal quantum number of the second orbital.
        a2 (float): Exponential decay of the second orbital.
        n3 (int): Principal quantum number of the third orbital.
        a3 (float): Exponential decay of the third orbital.
        n4 (int): Principal quantum number of the fourth orbital.
        a4 (float): Exponential decay of the fourth orbital.
        L (int): Laplace expansion order for Coulomb repulsion

    Returns:
        float: Evaluation of the radial integral.
    """
    fun = lambda r1, r2: r1**2 * r2**2 * r1**(n1 + n3 - 2) * r2**(n2 + n4 - 2) * np.exp(-(a1 + a3) * r1) * np.exp(-(a2 + a4) * r2) * np.min((r1, r2))**L / np.max((r1, r2))**(L + 1)

    return spi.dblquad(
        fun,
        0, 300,
        0, 300
    )

def symbolic_radial_integral(n1, a1, n2, a2, n3, a3, n4, a4, L):
    """Computes the radial integral of the Slater-type orbitals using symbolic integration.

    Args:
        n1 (int): Principal quantum number of the first orbital.
        a1 (float): Exponential decay of the first orbital.
        n2 (int): Principal quantum number of the second orbital.
        a2 (float): Exponential decay of the second orbital.
        n3 (int): Principal quantum number of the third orbital.
        a3 (float): Exponential decay of the third orbital.
        n4 (int): Principal quantum number of the fourth orbital.
        a4 (float): Exponential decay of the fourth orbital.
        L (int): Laplace expansion order for Coulomb repulsion

    Returns:
        float: Evaluation of the radial integral.
    """

    r1, r2 = sp.symbols('r1 r2')
    fun = r1**(n1 + n3 - 2) * r2**(n2 + n4 - 2) * sp.exp(-(a1 + a3) * r1) * sp.exp(-(a2 + a4) * r2) * sp.Min(r1, r2)**L / sp.Max(r1, r2)**(L + 1)

    return sp.integrate(
        fun * r1**2 * r2**2,
        (r1, 0, sp.oo),
        (r2, 0, sp.oo)
    )

def numerical_angular_integral(l1, m1, l2, m2, L, M):
    """Computes the angular integral of the Slater-type orbitals using numerical integration.

    Args:
        l1 (int): First azimuthal quantum number.
        m1 (int): First magnetic quantum number.
        l2 (int): Second azimuthal quantum number.
        m2 (int): Second magnetic quantum number.
        L (int): Radial order for Coulomb repulsion.
        M (int): Laplace expansion coefficient.

    Returns:
        float: Evaluation of the angular integral.
    """

    fun = lambda phi, theta: np.sin(theta) * np.conj(sps.sph_harm(m1, l1, phi, theta)) * sps.sph_harm(m2, l2, phi, theta) * sps.sph_harm(M, L, phi, theta)
    real_fun = lambda phi, theta: np.real(fun(phi, theta))
    imaginary_fun = lambda phi, theta: np.imag(fun(phi, theta))

    real_integral = spi.dblquad(
        real_fun,
        0, np.pi,
        0, 2 * np.pi
    )

    imaginary_integral = spi.dblquad(
        imaginary_fun,
        0, np.pi,
        0, 2 * np.pi
    )

    return (real_integral[0], imaginary_integral[0])

def symbolic_angular_integral(l1, m1, l2, m2, L, M):
    """Computes the angular integral of the Slater-type orbitals using symbolic integration.

    Args:
        l1 (int): First azimuthal quantum number.
        m1 (int): First magnetic quantum number.
        l2 (int): Second azimuthal quantum number.
        m2 (int): Second magnetic quantum number.
        L (int): Radial order for Coulomb repulsion.
        M (int): Laplace expansion coefficient.

    Returns:
        float: Evaluation of the angular integral.
    """

    phi, theta = sp.symbols('phi theta')

    Y1 = sh.Ynm(l1, m1, theta, phi).expand(func=True)
    Y2 = sh.Ynm(l2, m2, theta, phi).expand(func=True)
    Y3 = sh.Ynm(L, M, theta, phi).expand(func=True)

    integrand = sp.conjugate(Y1) * Y2 * Y3
    integral = sp.integrate(sp.sin(theta) * integrand, (theta, 0, sp.pi), (phi, 0, 2 * sp.pi))

    return integral

def compute_radial_integral(n1, l1, a1, n2, l2, a2, n3, l3, a3, n4, l4, a4, enable_numerical=False):
    print("========== Radial Integral ==========")
    print(f"n1 = {n1}, l1 = {l1}, a1 = {a1}")
    print(f"n2 = {n2}, l2 = {l2}, a2 = {a2}")
    print(f"n3 = {n3}, l3 = {l3}, a3 = {a3}")
    print(f"n4 = {n4}, l4 = {l4}, a4 = {a4}")

    L = 0
    Lmax = np.min((l1 + l3, l2 + l4))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        while L <= Lmax:
            print(f"Symbolic (L={L}): {symbolic_radial_integral(n1, a1, n2, a2, n3, a3, n4, a4, L).n():.15e}")
            if enable_numerical:
                print(f"Numerical (L={L}): {numerical_radial_integral(n1, a1, n2, a2, n3, a3, n4, a4, L).n():.15e}")
            L += 1

    print("========== End Radial Integral ==========")
    print("")

def compute_angular_integral(l1, m1, l2, m2, enable_numerical=False):
    print("========== Angular Integral ==========")
    print(f"l1 = {l1}, m1 = {m1}")
    print(f"l2 = {l2}, m2 = {m2}")

    Lmax = l1 + l2
    M = m1 - m2
    L = abs(M)

    if (abs(M) > Lmax):
        print("All angular integrals are zero.")
        Lmax = -1
    else:
        print(f"Angular integrals non zero <=> M = {M}")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        while L <= Lmax:
            symbolic_result = sp.N(symbolic_angular_integral(l1, m1, l2, m2, L, M))
            if symbolic_result != 0:
                print(f"Symbolic (L={L}): {symbolic_result:.15e}")
            else:
                print(f"Symbolic (L={L}): 0.0")

            if enable_numerical:
                print(f"Numerical: {numerical_angular_integral(l1, m1, l2, m2, L, M):.15e}")

            L += 1

    print("========== End Angular Integral ==========")
    print("")

if __name__ == "__main__":
    compute_radial_integral(2, 1, 1, 3, 1, 1, 2, 1, 1, 2, 1, 1)
    compute_radial_integral(2, 1, 3, 3, 1, 6, 3, 1, 9, 2, 1, 12)
    compute_radial_integral(4, 3, 3, 3, 2, 6, 5, 4, 9, 3, 2, 12)
    compute_radial_integral(2, 1, 50, 3, 1, 47, 2, 1, 50, 3, 1, 47)

    compute_angular_integral(0, 0, 0, 0)
    compute_angular_integral(2, 0, 0, 0)
    compute_angular_integral(2, 0, 1, 0)
    compute_angular_integral(3, -2, 4, 2)
    compute_angular_integral(3, -2, 4, 3)
    compute_angular_integral(3, -2, 3, 2)