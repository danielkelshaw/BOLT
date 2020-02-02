#ifndef UTILS_H
#define UTILS_H

#include "input.h"

namespace Utils {

double extrapolate(const std::vector<double> &vec, const std::vector<int> &normal, int order, int i, int j, int d = 0, int dims = 1);
double diracDelta(double dist);

}

#endif
