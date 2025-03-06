#pragma once

#include "Orbitals/orbital.hpp"

/**
 * @brief Slater type orbital (STO) representation

    It is the braket representation \f$ | R_n^\alpha l m \rangle \f$ where :
    \f[
        R_n^\alpha (r) = r^{n-1}  e^{-\alpha r}
    \f]
    and \f$ | l m \rangle \f$ are angular part of the wavefunction which are the
 spherical harmonics.
 **/
class SlaterPrimitive final : public Orbital {
  public:
    SlaterPrimitive() = default;
    ~SlaterPrimitive() = default;

    /**
     * @note The normalization constant only normalize the radial part of the
     * wavefunction as we expressed braket STO without normalization constant.
     * The angular part is normalized in the spherical harmonics.
     *
     * @throws std::invalid_argument if one or more parameters are invalid.
     *
     * @param n Principal quantum number (n > 0)
     * @param l Secondary quantum number (0 < l < n)
     * @param m Magnetic quantum number (-l <= m <= l)
     * @param alpha Exponential decay constant (alpha > 0)
     */
    SlaterPrimitive(int n, int l, int m, double alpha);

    /**
     * @note Updates the normalization constant if parameter is accepted.
     *
     * @throws std::invalid_argument if one or more parameters are invalid.
     *
     * @param alpha Exponential decay constant (alpha > 0)
     */
    void set_alpha(double alpha);

    /**
     * @note Updates the normalization constant if parameter is accepted.
     *
     * @throws std::invalid_argument if one or more parameters are invalid.
     *
     * @param n Principal quantum number (n > 0)
     */
    void set_n(int n);

    /**
     * @note Updates the normalization constant if parameter is accepted.
     *
     * @throws std::invalid_argument if one or more parameters are invalid.
     *
     * @param l Azimuthal quantum number (0 < l < n)
     */
    void set_l(int l);

    /**
     * @note Updates the normalization constant if parameter is accepted.
     *
     * @throws std::invalid_argument if one or more parameters are invalid.
     *
     * @param m Magnetic quantum number (-l <= m <= l)
     */
    void set_m(int m);

    inline int n() const { return m_n; }

    inline int l() const { return m_l; }

    inline int m() const { return m_m; }

    inline double alpha() const { return m_alpha; }

    friend std::ostream &operator<<(std::ostream &os,
                                    const SlaterPrimitive &orbital) {
        os << "(";
        os << orbital.n() << ", ";
        os << orbital.l() << ", ";
        os << orbital.m() << ") exp(- ";
        os << orbital.alpha() << " r)";
        return os;
    }

  private:
    void check_parameters(const int n, const int l, const int m,
                          const double alpha);
    void update_constant();

    int m_n, m_l, m_m;
    double m_alpha;
};

typedef SlaterPrimitive STO;

/**
 * @brief Computes the overlap integral between two slater primitives <o1|o2>
 *
 * @param o1 First slater primitive
 * @param o2 Second slater primitive
 * @return double <o1|o2>
 */
double overlap_integral(const SlaterPrimitive &, const SlaterPrimitive &,
                        const int n_offset = 0);

/**
 * @brief Computes the laplacian integral between two slater primitives
 *
 * @param o1 First slater primitive
 * @param o2 Second slater primitive
 * @return double <o1|nabla^2|o2>
 */
double laplacian_integral(const SlaterPrimitive &, const SlaterPrimitive &);

/**
 * @brief Computes the 1/|r-pos| integral between two slater primitives
 *
 * @param o1 First slater primitive
 * @param o2 Second slater primitive
 * @param pos Position of the considered nucleus
 * @return double <o1|1/|r-pos| |o2>
 */
double electron_nucleus_integral(const SlaterPrimitive &,
                                 const SlaterPrimitive &,
                                 const Eigen::Vector3d &);

/**
 * @brief Computes the electron-electron integral between four slater
 * primitives
 *
 * @param o1 First slater primitive
 * @param o2 Second slater primitive sharing coordinates of o1
 * @param o3 Third slater primitive
 * @param o4 Fourth slater primitive sharing coordinates of o3
 * @return double <o1 o3 | 1/|r1-r2| | o2 o4>
 */
double electron_electron_integral(const SlaterPrimitive &,
                                  const SlaterPrimitive &,
                                  const SlaterPrimitive &,
                                  const SlaterPrimitive &);