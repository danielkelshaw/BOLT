// Includes
#include "../include/Utils.h"

double Utils::extrapolate(const std::vector<double> &vec, const std::vector<int> &normal, int order, int i, int j, int d, int dims) {

	if (order == 0) {

		int i1 = i + normal[eX];
		int j1 = j + normal[eY];

		return vec[(i1 * Ny + j1) * dims + d];
	}
	else if (order == 1) {

		int i1 = i + normal[eX];
		int i2 = i + 2 * normal[eX];
		int j1 = j + normal[eY];
		int j2 = j + 2 * normal[eY];

		return 2.0 * vec[(i1 * Ny + j1) * dims + d] - vec[(i2 * Ny + j2) * dims + d];
	}
	else if (order == 2) {

		int i1 = i + normal[eX];
		int i2 = i + 2 * normal[eX];
		int i3 = i + 3 * normal[eX];
		int j1 = j + normal[eY];
		int j2 = j + 2 * normal[eY];
		int j3 = j + 3 * normal[eY];

		return 2.5 * vec[(i1 * Ny + j1) * dims + d] - 2.0 * vec[(i2 * Ny + j2) * dims + d] + 0.5 * vec[(i3 * Ny + j3) * dims + d];
	}
	else {
		std::cout << std::endl << std::endl << "Invalid Order" << std::endl << std::endl;
		exit(-1);
	}
}
