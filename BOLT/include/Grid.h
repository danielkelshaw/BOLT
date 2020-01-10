#ifndef GRID_H
#define GRID_H

// Includes
#include "input.h"

class GridClass {
public:
	GridClass();
	~GridClass() {};

public:

	// Grid parameters
	int t;											// Time step

	// Scaling parameters
	double Dx;										// Length scaling
	double Dt;										// Time scaling
	double Dm;										// Mass scaling
	double Drho;									// Density scaling

private:

	// Lattice parameters
	int nVels;										// Number of velocities
	double c_s;										// Speed of sound
	std::vector<double> w;							// Weightings
	std::vector<std::vector<int> > c;				// Direction vectors
	std::vector<int> opposite;						// Opposite vectors

	// Grid parameters
	double tau;										// Relaxation time
	double nu;										// Viscosity

	// Flattened kernel arrays
	std::vector<double> u;							// Velocity
	std::vector<double> u_n;						// Velocity (start of timestep)
	std::vector<double> rho;						// Density
	std::vector<double> rho_n;						// Density (start of timestep)
	std::vector<double> f;							// Populations
	std::vector<double> f_n;						// Populations (start of timestep)
	std::vector<eLatType> type;						// Lattice type matrix

	// Boundary conditions
	std::vector<int> BCVec;							// Site IDs to apply BCs
	std::vector<std::vector<double>> delU;			// Convective speed through boundary

	// Other
	std::vector<double> u_in;						// Inlet velocity profile

public:

	// LBM Methods
	void solver();									// Main solver

private:

	// LBM Methods
	void lbmKernel();												// Main LBM Kernel
	void streamCollide(int i, int j, int id);						// Stream/Collide: Pull algorithm

	double equilibrium(int id, int v);								// Equilibrium function
	double latticeForce(int id, int v);								// Discretise lattice force
	void macroscopic(int id);										// Calculate macroscopic quantities

	void applyBC(int i, int j, int id);								// Apply BCs
	void convectiveBC(int j, int id);								// Calculates f at boundary for eConvective
	void convectiveSpeed();											// Calculate convective speed
	void regularisedBC(int i, int j, int id,
			std::vector<int> &normalVector, eDirectionType normal);	// Regularised BC

	void initialiseGrid();											// Initialise the grid

	std::vector<int> getNormalVector(int i, int j, eDirectionType &normalDirection);

};

#endif
