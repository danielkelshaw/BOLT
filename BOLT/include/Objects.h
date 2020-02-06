#ifndef OBJECTS_H
#define OBJECTS_H

// Includes
#include "Grid.h"
#include "IBMBody.h"
#include "Utils.h"
#include "Overloads.h"

// Forward declaration
class GridClass;
class IBMBodyClass;

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

    GridClass *gridPtr;    // pointer to grid
    bool hasIBM;        // flag say if IBM bodies

// private members
private:

    std::vector<std::vector<size_t>> idxIBM; // vector containing body/node index for every node
    std::vector<IBMBodyClass> ibmBody;

    // Sub-iteration parameters
    int subIt;
    double relax;
    double subRes;
    double subNum;
    double subDen;

// public methods
public:

    // object routines
    void objectKernel();    // main kernel for objects

    void writeVTK();        // write IBM VTK file

// private methods
private:

    // IBM routines
    void ibmKernel();       // basic calculations

    // initialisation
    void initialiseObjects();
    void computeEpsilon();

    // input / output
    void readGeometry();
    void resetPointers();
};

#endif
