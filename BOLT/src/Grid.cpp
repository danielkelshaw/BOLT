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
	type.resize(Nx * Ny, eFluid);

	u_in.resize(Ny * dims, 0.0);

	initialiseGrid();
}

void GridClass::initialiseGrid() {
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {

			int id = i * Ny + j;

			// Left wall
			if (i == 0) {
				type[id] = eVelocity;
			}

			// Right wall
			if (i == Nx - 1) {
				type[id] = eConvective;
			}

			// Bottom wall
			if (j == 0) {
				type[id] = eVelocity;
			}

			// Top wall
			if (j == Ny - 1) {
				type[id] = eVelocity;
			}

			if (type[id] != eFluid) {
				BCVec.push_back(id);
			}
		}

		// Set inlet velocity
		for (int j = 0; j < Ny; j++) {
			u_in[j * dims + eX] = ux0_p * Dt / Dx;
			u_in[j * dims + eY] = uy0_p * Dt / Dx;
		}
	}

	// Set initial velocity
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {

			int id = i * Ny + j;

			u[id * dims + eX] = ux0_p * Dt / Dx;
			u[id * dims + eY] = uy0_p * Dt / Dx;

			if (type[id] == eWall) {
				u[id * dims + eX] = 0.0;
				u[id * dims + eY] = 0.0;
			}
		}
	}

	u_n = u;
	rho_n = rho;

	// Set f values to equilibrium
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {

			int id = i * Ny + j;

			for (int v = 0; v < nVels; v++) {
				f[id * nVels + v] = equilibrium(id, v);
			}
		}
	}

	f_n = f;
}

double GridClass::equilibrium(int id, int v) {
	return 0.0;
}
