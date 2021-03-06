// Includes
#include "../include/IBMBody.h"


IBMBodyClass::IBMBodyClass(ObjectsClass *objects, int bodyID, const std::vector<double> &pos, double radius) {

    objectPtr = objects;
    ID = bodyID;
    bodyType = eCircle;

    // count number of nodes
    int numNodes = static_cast<int>(std::floor(2.0 * M_PI * radius / objectPtr->gridPtr->Dx));

    // at the moment this position doesnt do anything - will be passed to a node class
    std::vector<double> position(dims);

    for (int i = 1; i < numNodes; i++) {

        // get position of body
        position[eX] = pos[eX] + radius * std::cos(i * 2.0 * M_PI / numNodes);
        position[eY] = pos[eY] + radius * std::sin(i * 2.0 * M_PI / numNodes);

        // construct nodes
        nodes.emplace_back(this, (i - 1), position);
    }
}
