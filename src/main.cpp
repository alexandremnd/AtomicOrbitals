#include "Atom/atom.hpp"
#include "BasisSet/slater_primitive.hpp"
#include "Eigen/Dense"
#include <algorithm>
#include <iostream>
#include <memory>
#include <vector>

int main() {
    auto He = Atom<SlaterPrimitive>(2, Eigen::Vector3d(0, 0, 0));
    He.add_orbital(2, 0, 0, 1.3);


    // RestrictedHartreeFock<SlaterPrimitive> hf(atom);

    return 0;
}