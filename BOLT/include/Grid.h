#ifndef GRID_H
#define GRID_H

// Includes
#include "input.h"
#include "Objects.h"

// Forward declaration
class ObjectsClass;

class GridClass {

    // Friend classes
    friend class ObjectsClass;
    friend class IBMNodeClass;

public:

    GridClass();              // Default constructor
    ~GridClass() {};          // Default destructor

public:

    // Grid parameters
    int t;                    // Time step

    // Object Pointer
    ObjectsClass *objectPtr;

    // Scaling parameters
    double Dx;                // Length scaling
    double Dt;                // Time scaling
    double Dm;                // Mass scaling
    double Drho;              // Density scaling

private:

    // Lattice parameters
    int nVels;                              // Number of velocities
    double c_s;                             // Speed of sound
    std::vector<double> w;                  // Weightings
    std::vector<std::vector<int> > c;       // Direction vectors
    std::vector<int> opposite;              // Opposite vectors

    // Grid parameters
    double tau;         // Relaxation time
    double nu;          // Viscosity

    // Flattened kernel arrays
    std::vector<double> u;                  // Velocity
    std::vector<double> u_n;                // Velocity (start of timestep)
    std::vector<double> rho;                // Density
    std::vector<double> rho_n;              // Density (start of timestep)

    std::vector<double> omegaLocal;         // Local omega value
    std::vector<double> force_xyz;          // Cartestian force (pressure and gravity)
    std::vector<double> force_ibm;          // Cartestian force (IBM)

    std::vector<double> f;                  // Populations
    std::vector<double> f_n;                // Populations (start of timestep)
    std::vector<eLatType> type;             // Lattice type matrix

    // Boundary conditions
    std::vector<int> BCVec;                 // Site IDs to apply BCs
    std::vector<std::vector<double>> delU;  // Convective speed through boundary

    // Other
    std::vector<double> u_in;               // Inlet velocity profile

    // Timings
    std::chrono::time_point<std::chrono::steady_clock> startTime;   // Time when clock is called
    std::chrono::duration<double> loopTime;                         // Average loop time

public:

    // LBM Methods
    void solver();      // Main solver
    void writeInfo();   // Write info to screen
    void writeVTK();    // Write info to .VTK file

private:

    // LBM Methods
    void lbmKernel();                              // Main LBM Kernel
    void streamCollide(int i, int j, int id);      // Stream/Collide: Pull algorithm

    double equilibrium(int id, int v);             // Equilibrium function
    double latticeForce(int id, int v);            // Discretise lattice force
    void macroscopic(int id);                      // Calculate macroscopic quantities

    // Boundary Condition Methods
    void applyBC(int i, int j, int id);            // Apply BCs
    void convectiveBC(int j, int id);              // Calculates f at boundary for eConvective
    void convectiveSpeed();                        // Calculate convective speed
    void regularisedBC(int i, int j, int id,
            std::vector<int> &normalVector, eDirectionType normal);   // Regularised BC

    void initialiseGrid();                         // Initialise the grid

    // Helper Methods
    std::vector<int> getNormalVector(int i, int j, eDirectionType &normalDirection);

};

#endif
