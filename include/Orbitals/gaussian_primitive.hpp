#pragma once

#include "Eigen/Dense"
#include "Orbitals/orbital.hpp"

/**
 * @brief Gaussian primitive orbital.
 * Orbital defined as:
 * \f[
 *      \phi(x, y, z) = N \cdot (x-x_A)^i (y-y_A)^j (z-z_A)^k \cdot e^{-\alpha
 *      ||\mathbf{r} - \mathbf{R_A}||^2}
 * \f]
 * where \f$N\f$ is the normalization constant, \f$\alpha\f$ is the decay rate,
 * \f$i, j, k\f$ are the exponents in the x, y, z directions respectively, and
 * \f$\mathbf{R_A}\f$ is the position of the nucleus.
 *
 */
class GaussianPrimitive : public Orbital {
  public:
    GaussianPrimitive(){};

    /**
     * @brief Construct a new Gaussian Primitive object
     *
     * @param x_exponent Power for the \f$(x-x_A)^i \f$ factor
     * @param y_exponent Power for the \f$(y-y_A)^j \f$ factor
     * @param z_exponent Power for the \f$(z-z_A)^k \f$ factor
     * @param alpha Decay rate of the exponential \f$ e^{-\alpha
     *      ||\mathbf{r} - \mathbf{R_A}||^2} \f$
     * @param position Coordinates of the \f$ \mathbf{R_A} = (x_A, y_A, z_A) \f$
     * vector which is the center of the orbital (often being the nucleus
     * position)
     */
    GaussianPrimitive(int x_exponent, int y_exponent, int z_exponent,
                      double alpha, Eigen::Vector3d position = {0, 0, 0}) {
        m_position = position;
        set_x_exponent(x_exponent);
        set_y_exponent(y_exponent);
        set_z_exponent(z_exponent);
        set_alpha(alpha);
        update_normalisation();
    }

    /**
     * @brief Evaluates the orbital at a given position.
     *
     * @param x x-coordinate
     * @param y y-coordinate
     * @param z z-coordinate
     * @return double \f$ \Phi(x, y, z) \f$
     */
    double evaluate(double x, double y, double z) const;

    void update_normalisation();

    // Setters & Getters
    /**
     * @throw std::invalid_argument if alpha <= 0
     */
    void set_alpha(double alpha);

    /**
     * @throw std::invalid_argument if x_exponent < 0
     */
    void set_x_exponent(int x_exponent);

    /**
     * @throw std::invalid_argument if y_exponent < 0
     */
    void set_y_exponent(int y_exponent);

    /**
     * @throw std::invalid_argument if z_exponent < 0
     */
    void set_z_exponent(int z_exponent);

    int x_exponent() const { return m_x_exponent; }
    int y_exponent() const { return m_y_exponent; }
    int z_exponent() const { return m_z_exponent; }
    double alpha() const { return m_alpha; }

    friend GaussianPrimitive operator*(const GaussianPrimitive &orbital1,
                                       const GaussianPrimitive &orbital2);

    GaussianPrimitive &operator=(const GaussianPrimitive &rhs);

    friend std::ostream &operator<<(std::ostream &os,
                                    const GaussianPrimitive &orbital) {
        os << "(";
        os << orbital.m_x_exponent << ", ";
        os << orbital.m_y_exponent << ", ";
        os << orbital.m_z_exponent << ") exp(- ";
        os << orbital.m_alpha << " ||r - R||^2)";
        return os;
    }

  private:
    int m_x_exponent;
    int m_y_exponent;
    int m_z_exponent;
    double m_alpha;
};

typedef GaussianPrimitive GTO;

/**
 * @brief Computes the overlap integral between two Gaussian primitives <o1|o2>
 *
 * @param o1 First gaussian primitive
 * @param o2 Second gaussian primitive
 * @return double <o1|o2>
 */
double overlap_integral(const GaussianPrimitive &o1,
                        const GaussianPrimitive &o2);

/**
 * @brief Computes the laplacian integral between two Gaussian primitives
 *
 * @param o1 First gaussian primitive
 * @param o2 Second gaussian primitive
 * @return double <o1|nabla^2|o2>
 */
double laplacian_integral(const GaussianPrimitive &o1,
                          const GaussianPrimitive &o2);

/**
 * @brief Computes the 1/|r-pos| integral between two Gaussian primitives
 *
 * @param o1 First gaussian primitive
 * @param o2 Second gaussian primitive
 * @param pos Position of the considered nucleus
 * @return double <o1|1/|r-pos| |o2>
 */
double electron_nucleus_integral(const GaussianPrimitive &o1,
                                 const GaussianPrimitive &o2,
                                 const Eigen::Vector3d &pos);

/**
 * @brief Computes the electron-electron integral between four Gaussian
 * primitives
 *
 * @param o1 First gaussian primitive
 * @param o2 Second gaussian primitive sharing coordinates of o1
 * @param o3 Third gaussian primitive
 * @param o4 Fourth gaussian primitive sharing coordinates of o3
 * @return double <o1 o3 | 1/|r1-r2| | o2 o4>
 */
double electron_electron_integral(const GaussianPrimitive &o1,
                                  const GaussianPrimitive &o2,
                                  const GaussianPrimitive &o3,
                                  const GaussianPrimitive &o4);