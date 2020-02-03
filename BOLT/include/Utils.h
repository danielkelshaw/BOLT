#ifndef UTILS_H
#define UTILS_H

#include "input.h"
#include <Accelerate/Accelerate.h>

namespace Utils {

double extrapolate(const std::vector<double> &vec, const std::vector<int> &normal, int order, int i, int j, int d = 0, int dims = 1);
double diracDelta(double dist);

std::vector<double> solveLAPACK(std::vector<double> A, std::vector<double> b, int BC = 0);

}

#endif
