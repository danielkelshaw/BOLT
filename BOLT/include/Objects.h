#ifndef OBJECTS_H
#define OBJECTS_H

// Includes
#include "Grid.h"

class ObjectsClass {

// constructors / destructors
public:

    // default constructor / destructor
    ObjectsClass();
    ~ObjectsClass() {};

    // custom constructor
    ObjectsClass(GridClass &grid);
    
// public members
public:

    GridClass *gPtr;    // pointer to grid
    bool hasIBM;        // flag say if IBM bodies

// private members
private:

    std::vector<std::vector<size_t>> idxIBM; // vector containing body/node index for every node

// public methods
public:

    // object routines
    void objectKernel();    // main kernel for objects

// private methods
private:

    // IBM routines
    void ibmKernel();       // basic calculations

    // initialisation
    void initialiseObjects();

    // input / output
    void readGeometry();
};

#endif