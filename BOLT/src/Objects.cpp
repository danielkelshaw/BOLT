// Includes
#include "../include/Objects.h"

void ObjectsClass::objectKernel() {

    // At the moment only considering IBM.
    ibmKernel();
}

void ObjectsClass::ibmKernel() {

    fill(gridPtr->force_ibm.begin(), gridPtr->force_ibm.begin(), 0.0);

    for (size_t i = 0; i < idxIBM.size(); i++) {

        size_t ib = idxIBM[i][0];
        size_t n = idxIBM[i][1];

        ibmBody[ib].node[n].interpolate();
        ibmBody[ib].node[n].forceCalc();
        ibmBody[ib].node[n].spread();
    }

    for (size_t i = 0; i < idxIBM.size(); i++) {

        size_t ib = idxIBM[i][0];
        size_t n = idxIBM[i][1];

        ibmBody[ib].node[n].updateMacroscopic();
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
        else {
            std::cout << "Only CIRCLE is supported at the moment..." << std::endl;
            exit(-1);
        }

        bodyCase = "NONE";
    }

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

    std::vector<IBMBodyClass> *ibmBodyPtr = &ibmBody;

    double Dx = gridPtr->Dx;

    // Loop through each body
    for (size_t b = 0; b < ibmBodyPtr->size(); b++) {

        if (gridPtr->t == 0) {

            size_t dim = (*ibmBodyPtr)[b].nodes.size();

            std::vector<double> A(dim * dim, 0.0);

            for (size_t i = 0; i < dim; i++) {
                for (size_t j = 0; j < dim; j++) {

                    for (size_t s = 0; s < (*ibmBodyPtr)[b].nodes[i].supps.size(); s++) {

                        double diracVal_i = (*ibmBodyPtr)[b].nodes[i].supps[s].diracVal;

                        double distX = fabs((*ibmBodyPtr)[b].nodes[j].pos[eX] / Dx - (*ibmBodyPtr)[b].nodes[i].supps[s].idx);
                        double distY = fabs((*ibmBodyPtr)[b].nodes[j].pos[eY] / Dx - (*ibmBodyPtr)[b].nodes[i].supps[s].jdx);
                        double diracVal_j = Utils::diracDelta(distX) * Utils::diracDelta(distY);

                        A[i * dims + j] += diracVal_i * diracVal_j;

                    }

                    A[i * dims + j] *= (*ibmBodyPtr)[b].nodes[j].ds;
                }
            }

            std::vector<double> b_vec(dim, 1.0);
            std::vector<double> epsilon = Utils::solveLAPACK(A, b_vec);

            for (size_t i = 0; i < dim; i++) {
                (*ibmBodyPtr)[b].nodes[i].epsilon = epsilon[i];
            }
        }
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
