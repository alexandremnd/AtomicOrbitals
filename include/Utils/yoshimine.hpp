#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>

/**
 * @brief Container using Yoshimine sorting.
 *
 * We use (ij | kl) as a shorthand for <ik | w | jl>.
 * For electron electron integrals denoted (ij | kl), we have an 8-fold
 * symmetry: (ij | kl) = (ji | kl) = (ij | lk) = (ji | lk) = (kl | ij) = (lk |
 * ij) = (kl | ji) = (lk | ji)
 *
 * For an unordered set of indices (i, j, k, l), index(i, j, k, l) = index(j, i,
 * k, l) = ...
 *
 * For a basis set of size N, Yoshimine container requires N(N+1)(N² + N + 2)/8
 * elements instead of N⁴ elements.
 */
template <typename T> class Yoshimine {
  public:
    /**
     * @param size Size of the Yoshimine container
     */
    Yoshimine(int size) {
        size = (size * (size + 1) * (size * size + size + 2)) / 8;
        m_yoshimine.resize(size);
    }

    Yoshimine() = default;

    /**
     * @brief Returns the corresponding element (ab| cd) in the Yoshimine
     * container. Used for reducing memory usage
     *
     * @param a First electron index
     * @param b Second electron index sharing the same coordinates as a
     * @param c Third electron index
     * @param d Fourth electron index sharing the same coordinates as c
     * @return T
     */
    inline T operator()(int a, int b, int c, int d) const {
        int abcd = index(a, b, c, d);

        if (abcd >= m_yoshimine.size()) {
            throw std::out_of_range("Yoshimine index out of range");
        }

        return m_yoshimine[abcd];
    }

    /**
     * @brief Returns the corresponding element (ab| cd) in the Yoshimine
     * container. Used for reducing memory usage
     *
     * @param a First electron index
     * @param b Second electron index sharing the same coordinates as a
     * @param c Third electron index
     * @param d Fourth electron index sharing the same coordinates as c
     * @return T
     */
    inline T &operator()(int a, int b, int c, int d) {
        int abcd = index(a, b, c, d);

        if (abcd >= m_yoshimine.size()) {
            throw std::out_of_range("Yoshimine index out of range");
        }

        return m_yoshimine[abcd];
    }

    void print_content() const {
        for (int i = 0; i < m_yoshimine.size(); i++) {
            std::cout << i << ": " << m_yoshimine[i] << std::endl;
        }
    }

    inline int size() const { return m_yoshimine.size(); }

  private:
    inline int index(int a, int b, int c, int d) const {
        int ab, cd, abcd;

        if (a > b) {
            ab = (a * (a + 1)) / 2 + b;
        } else {
            ab = (b * (b + 1)) / 2 + a;
        }

        if (c > d) {
            cd = (c * (c + 1)) / 2 + d;
        } else {
            cd = (d * (d + 1)) / 2 + c;
        }

        if (ab > cd) {
            abcd = (ab * (ab + 1)) / 2 + cd;
        } else {
            abcd = (cd * (cd + 1)) / 2 + ab;
        }

        return abcd;
    }

    std::vector<T> m_yoshimine;
};