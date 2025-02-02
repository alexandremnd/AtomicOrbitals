#include "BasisSet/gaussian_primitive.hpp"
#include <stdexcept>

double GaussianPrimitive::evaluate(double x, double y, double z) const {
    return 0.0;
}

void GaussianPrimitive::set_x_exponent(int x_exponent) {
    return;
}

void GaussianPrimitive::set_y_exponent(int y_exponent) {
    return;
}

void GaussianPrimitive::set_z_exponent(int z_exponent) {
    return;
}


// ===============================================================================================
// ============================ Matrix element for Gaussian Primitive ============================
// ===============================================================================================

double overlap_integral(const GaussianPrimitive& orbital1, const GaussianPrimitive& orbital2) {
    throw std::logic_error("overlap_integral not implemented for GaussianPrimitive");
}

double laplacian_integral(const GaussianPrimitive& orbital1, const GaussianPrimitive& orbital2) {
    throw std::logic_error("overlap_integral not implemented for GaussianPrimitive");
}

double electron_nucleus_integral(const GaussianPrimitive& orbital1, const GaussianPrimitive& orbital2, const Eigen::Vector3d& nucleus_position) {
    throw std::logic_error("overlap_integral not implemented for GaussianPrimitive");
}

double electron_electron_integral(const GaussianPrimitive &orbital1, const GaussianPrimitive &orbital2, const GaussianPrimitive &orbital3, const GaussianPrimitive &orbital4)
{
    throw std::logic_error("overlap_integral not implemented for GaussianPrimitive");
}