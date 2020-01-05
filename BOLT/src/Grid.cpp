// Includes
#include "../include/Grid.h"

GridClass::GridClass() {

	// Output message
	std::cout << "Initialising GridClass." << std::endl;

	// Default parameters
	double rho0 = 1.0;

	// D2Q9 parameters
	nVels = 9;
	c_s = 1.0 / std::sqrt(3.0);
	w = {4.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};
	c = {{0, 0}, {1, 0}, {-1, 0}, {0, 1}, {0, -1}, {1, 1}, {-1, -1}, {1, -1}, {-1, 1}};
	opposite = {0, 2, 1, 4, 3, 6, 5, 8, 7};

	// Initialise parameters
	t = 0;
	tau = 1.0 / omega;
	nu = (tau - 0.5) * SQ(c_s);
	Dx = height_p / (Ny - 1);
	Dt = (nu / nu_p) * SQ(Dx);
	Dm = (rho_p / rho0) * CU(Dx);
	Drho = (rho_p / rho0);

	// Set array size and initialise
	u.resize(Nx * Ny * dims, 0.0);
	u_n.resize(Nx * Ny * dims, 0.0);
	rho.resize(Nx * Ny, rho0);
	rho_n.resize(Nx * Ny, rho0);
	f.resize(Nx * Ny * nVels, 0.0);
	f_n.resize(Nx * Ny * nVels, 0.0);

	initialiseGrid();
}

void GridClass::initialiseGrid() {
}