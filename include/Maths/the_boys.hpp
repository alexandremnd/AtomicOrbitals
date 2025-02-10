#pragma once

#include <boost/math/special_functions/gamma.hpp>
#include <cmath>

double TheBoys_Function(int &n, double &x);

std::vector<double> TheBoys_Recurrence(int &n, double &x);