import sympy as sp
import numpy as np
import symbolic_slater_braket as sb

def symbolic_electron_electron_integral(n1, l1, m1, a1,
                                        n2, l2, m2, a2,
                                        n3, l3, m3, a3,
                                        n4, l4, m4, a4):
    Lmax = np.max((l1 + l3, l2 + l4))

    result = 0;
    for L in range(Lmax + 1):
        m_result = 0
        for M in range(-L, L + 1):
            print(f"Calculating integral for L = {L}, M = {M}")
            angular_1 = sb.symbolic_angular_integral(l1, m1, l3, m3, L, M)
            angular_2 = sb.symbolic_angular_integral(l2, m2, l4, m4, L, -M)
            m_result += angular_1 * angular_2

            print(f"Angular 1: {angular_1.n()} (L = {L}, M = {M})")
            print(f"Angular 2: {angular_2.n()} (L = {L}, M = {M})")

        radial = sb.symbolic_radial_integral(n1, a1, n2, a2, n3, a3, n4, a4, L)
        m_result *= radial
        m_result *= 4 * np.pi / (2 * L + 1)
        result += m_result

    return result.n()

if __name__ == "__main__":
    print("[Warning] This computation may take some time (near 3 minutes in my case).")
    print(symbolic_electron_electron_integral(2, 1, 0, 2,
                                              3, 2, 0, 2,
                                              4, 2, -1, 2,
                                              2, 1, -1, 2))