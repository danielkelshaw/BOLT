// Includes
#include "../include/Objects.h"

void ObjectsClass::objectKernel() {

    // At the moment only considering IBM.
    ibmKernel();
}

void ObjectsClass::ibmKernel() {

    fill(gridPtr->force_ibm.begin(), gridPtr->force_ibm.end(), 0.0);

    for (size_t i = 0; i < idxIBM.size(); i++) {

        size_t ib = idxIBM[i][0];
        size_t n = idxIBM[i][1];

        ibmBody[ib].nodes[n].interpolate();
        ibmBody[ib].nodes[n].forceCalc();
        ibmBody[ib].nodes[n].spread();
    }

    for (size_t i = 0; i < idxIBM.size(); i++) {

        size_t ib = idxIBM[i][0];
        size_t n = idxIBM[i][1];

        ibmBody[ib].nodes[n].updateMacroscopic();
    }
}

void ObjectsClass::resetPointers() {

    for (size_t ib = 0; ib < ibmBody.size(); ib++) {
        for (size_t n = 0; n < ibmBody[ib].nodes.size(); n++) {
            ibmBody[ib].nodes[n].ibmPtr = &ibmBody[ib];
        }
    }
}

void ObjectsClass::readGeometry() {

    // open config file
    std::ifstream file;
    file.open("input/geometry.config", std::ios::in);

    // handle failure to open file
    if (!file.is_open()) {
        std::cout << "Error opening geometry config file..." << std::endl;
        exit(-1);
    }

    // skip comments in config file
    std::streamoff fileOffset;
    std::string line;
    file.seekg(std::ios::beg);

    do {
        fileOffset = file.tellg();
        std::getline(file, line);

    } while (line[0] == '#' && !file.eof());

    // move cursor over comments
    file.seekg(fileOffset, std::ios::beg);

    std::string bodyCase;

    // read in geometries
    int bodyID = 0;
    while (file) {

        // get type of body
        file >> bodyCase;

        int number;
        file >> number;

        std::vector<double> start(dims);
        file >> start[eX];
        file >> start[eY];

        std::vector<double> space(dims);
        file >> space[eX];
        file >> space[eY];

        if (bodyCase == "CIRCLE") {

            double radius;
            file >> radius;

            for (int i = 0; i < number; i++) {

                std::vector<double> pos = {start[eX] + i * space[eX], start[eY] + i * space[eY]};
                ibmBody.emplace_back(this, bodyID, pos, radius);
                bodyID++;
            }

        }

        bodyCase = "NONE";
    }

    resetPointers();

    if (ibmBody.size() > 0) {
        hasIBM = true;
    }

    if (hasIBM == true) {
        std::cout << "found " << ibmBody.size() << " object(s)" << std::endl;
    }
    else if (hasIBM == false) {
        std::cout << "no objects found" << std::endl;
    }

}

void ObjectsClass::initialiseObjects() {
    
    // Loop through all bodies
    for (size_t b = 0; b < ibmBody.size(); b++) {

        // Loop through all nodes
        for (size_t n = 0; n < ibmBody[b].nodes.size(); n++) {

            idxIBM.push_back({b, n});

            ibmBody[b].nodes[n].findSupport();
            ibmBody[b].nodes[n].computeDs();

        }
    }

    computeEpsilon();
}

void ObjectsClass::computeEpsilon() {

    double Dx = gridPtr->Dx;

    // Loop through each body
    for (size_t ib = 0; ib < ibmBody.size(); ib++) {

        if (gridPtr->t == 0) {

            size_t dim = ibmBody[ib].nodes.size();

            std::vector<double> A(dim * dim, 0.0);

            // Loop through all nodes
            for (size_t i = 0; i < dim; i++) {
                for (size_t j = 0; j < dim; j++) {

                    // Now loop through all support markers for node i
                    for (size_t s = 0; s < ibmBody[ib].nodes[i].supps.size(); s++) {

                        // Dirac delta value of support marker for node i and support s
                        double diracVal_i = ibmBody[ib].nodes[i].supps[s].diracVal;

                        double distX = fabs(ibmBody[ib].nodes[j].pos[eX] / Dx - ibmBody[ib].nodes[i].supps[s].idx);
                        double distY = fabs(ibmBody[ib].nodes[j].pos[eY] / Dx - ibmBody[ib].nodes[i].supps[s].jdx);
                        double diracVal_j = Utils::diracDelta(distX) * Utils::diracDelta(distY);

                        // Add to A matrix
                        A[i * dim + j] += diracVal_i * diracVal_j;
                    }

                    // Mulitply by volume
                    A[i * dim + j] *= ibmBody[ib].nodes[j].ds;
                }
            }

            std::vector<double> b_vec(dim, 1.0);
            std::vector<double> epsilon = Utils::solveLAPACK(A, b_vec);

            for (size_t i = 0; i < dim; i++) {
                ibmBody[ib].nodes[i].epsilon = epsilon[i];
            }
        }
    }
}

void ObjectsClass::writeVTK() {

    if (hasIBM == true) {

        std::ofstream output;
        output.precision(PRECISION);
        output.open("Results/VTK/IBM." + std::to_string(gridPtr->t) + ".vtk");

        // Write VTK header
        output << "# vtk DataFile Version 3.0\n";
        output << "IBM\n";
        output << "ASCII\n";
        output << "DATASET POLYDATA\n";

        size_t totalNodes = 0;
        size_t totalElements = 0;

        for (size_t ib = 0; ib < ibmBody.size(); ib++) {

            totalNodes += ibmBody[ib].nodes.size();

            if (ibmBody[ib].bodyType == eCircle) {
                totalElements += ibmBody[ib].nodes.size();
            }
        }

        output << "POINTS " << totalNodes << " double\n";

        for (size_t ib = 0; ib < ibmBody.size(); ib++) {
            for (size_t n = 0; n < ibmBody[ib].nodes.size(); n++) {
                output << ibmBody[ib].nodes[n].pos[eX] << " " << ibmBody[ib].nodes[n].pos[eY] << " 1\n";
            }
        }

        output << "LINES " << totalElements << " " << 3 * totalElements << "\n";

        size_t bodyOffset = 0;
        for (size_t ib = 0; ib < ibmBody.size(); ib++) {
            size_t nElements = 0;

            if (ibmBody[ib].bodyType == eCircle) {
                nElements = ibmBody[ib].nodes.size();
            }

            for (size_t el = 0; el < nElements; el++) {
                output << "2 " << bodyOffset + (el % ibmBody[ib].nodes.size())  << " " << bodyOffset + ((el + 1) % ibmBody[ib].nodes.size()) << "\n";
            }
            bodyOffset += ibmBody[ib].nodes.size();
        }

        output.close();
    }
}

ObjectsClass::ObjectsClass(GridClass &grid) {

    std::cout << "Initialising Objects..." << std::endl;

    // Set pointers
    gridPtr = &grid;
    gridPtr->objectPtr = this;

    // Initial sub-iteration values
    subIt = 0;
    subRes = 0.0;
    relax = 1.0;
    subNum = 0.0;
    subDen = 0.0;

    hasIBM = false;

    // Read in geometry
    readGeometry();

    // Initialise objects
    initialiseObjects();
}
