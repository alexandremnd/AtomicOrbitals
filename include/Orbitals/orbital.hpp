
#include "Eigen/Dense"

/**
 * @brief Base class for orbitals.
 * @note Basis orbitals must inherit from this class for the Hartree-Fock
 * method.
 */
class Orbital {
  public:
    ~Orbital() = default;

    double constant() const { return m_constant; };

    Eigen::Vector3d position() const { return m_position; };

    /**
     * @brief Set the orbital center position (nucleus position).
     * @note May be overriden by derived class.
     *
     * @param position New position
     */
    virtual void set_position(Eigen::Vector3d position) {
        m_position = position;
    };

  protected:
    double m_constant = 1;
    Eigen::Vector3d m_position = Eigen::Vector3d{0, 0, 0};
};