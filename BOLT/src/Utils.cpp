// Includes
#include "../include/Utils.h"

double Utils::extrapolate(const std::vector<double> &vec, const std::vector<int> &normal, int order, int i, int j, int d, int dims) {

    // Implements basic methods for extrapolation

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

double Utils::diracDelta(double dist) {

    double absDist = fabs(dist);

    if (absDist > 1.5) {
        return 0.0;
    }
    else if (absDist > 0.5) {
        return (5.0 - 3.0 * absDist - sqrt(-3.0 * SQ(1.0 - absDist) + 1.0)) / 6.0;
    }
    else {
        return (1.0 + sqrt(1.0 - 3.0 * SQ(absDist))) / 3.0;
    }
}

std::vector<double> Utils::solveLAPACK(std::vector<double> A, std::vector<double> b, int BC) {

    // Set up the correct values
    char trans = 'T';
    
    int dim = static_cast<int>(sqrt(static_cast<int>(A.size())));

    int row = dim - BC;
    int col = dim - BC;

    int offset = BC * dim + BC;
    int nrhs = 1;

    int LDA = dim;
    int LDB = dim;

    int info;
    std::vector<int> ipiv(row, 0);

    // Factorise and solve
    dgetrf_(&row, &col, A.data() + offset, &LDA, ipiv.data(), &info);
    dgetrs_(&trans, &row, &nrhs, A.data() + offset, &LDA, ipiv.data(), b.data() + BC, &LDB, &info);

    // Set return values not included to zero
    fill(b.begin(), b.begin() + BC, 0.0);

    // Return RHS
    return b;
}
