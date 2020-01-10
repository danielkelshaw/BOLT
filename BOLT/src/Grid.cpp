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

	// Assuming right wall is eConvective for now...
	delU.resize(Ny, std::vector<double>(dims, 0.0));

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

	double C1 = SQ(u_n[id * dims + eX]) + SQ(u_n[id * dims + eY]);
	double C2 = u_n[id * dims + eX] * c[v][eX] + u_n[id * dims + eY] * c[v][eY];

	return rho_n[id] * w[v] * (1.0 + 3.0 * C2 + 4.5 * SQ(C2) - 1.5 * C1);
}

void GridClass::streamCollide(int i, int j, int id) {

	for (int v = 0; v < nVels; v++) {
		int src_id = ((i - c[v][eX] + Nx) % Nx) * Ny + ((j - c[v][eY] + Ny) % Ny);
		f[id * nVels + v] = f_n[src_id * nVels + v] + omega * (equilibrium(src_id, v) - f_n[src_id * nVels + v]) + latticeForce(src_id, v);
	}
}

double GridClass::latticeForce(int id, int v) {
	// At the moment latticeForce is neglected for simplicity.
	return 0.0;
}

void GridClass::lbmKernel() {
	
	convectiveSpeed();

	f_n.swap(f);
	u_n.swap(u);
	rho_n.swap(rho);

	// Iterate through all points
	for (int i = 0; i < Nx; i++) {
		for (int j = 0; j < Ny; j++) {

			int id = i * Nx + j;

			streamCollide(i, j, id);

			if (type[id] == eFluid) {
				macroscopic(id);
			}
		}
	}

	// Apply BCs to all points
	for (int b = 0; b < BCVec.size(); b++) {

		int id = BCVec[b];
		int i = id / Ny;
		int j = id - (i * Ny);

		applyBC(i, j, id);
		macroscopic(id);
	}
}

void GridClass::convectiveSpeed() {

	// Set u_out to zero
	double u_out = 0.0;

	// Iterate through outlet cells
	for (int j = 0; j < Ny; j++) {
		u_out = u[((Nx - 1) * Ny + j) * dims + eX];
	}

	// Average result
	u_out /= static_cast<double>(Ny);

	// Set delU
	for (int j = 0; j < Ny; j++) {
		delU[j][eX] = (-u_out / 2.0) * (3.0 * u[((Nx - 1) * Ny + j) * dims + eX] - 4.0 * u[((Nx - 2) * Ny + j) * dims + eX] + u[((Nx - 3) * Ny + j) * dims + eX]);
		delU[j][eY] = (-u_out / 2.0) * (3.0 * u[((Nx - 1) * Ny + j) * dims + eY] - 4.0 * u[((Nx - 2) * Ny + j) * dims + eY] + u[((Nx - 3) * Ny + j) * dims + eY]);
	}
}

void GridClass::convectiveBC(int j, int id) {

	// Setting the f values assuming right wall is eConvective
	f[id * nVels + 2] = f_n[id * nVels + 2] + 3.0 * w[2] * (delU[j][eX] * c[2][eX] + delU[j][eY] * c[2][eY]);
	f[id * nVels + 6] = f_n[id * nVels + 6] + 3.0 * w[6] * (delU[j][eX] * c[6][eX] + delU[j][eY] * c[6][eY]);
	f[id * nVels + 8] = f_n[id * nVels + 8] + 3.0 * w[8] * (delU[j][eX] * c[8][eX] + delU[j][eY] * c[8][eY]);
}

void GridClass::macroscopic(int id) {

	// Reset values
	rho[id] = 0.0;
	u[id * dims + eX] = 0.0;
	u[id * dims + eY] = 0.0;

	// Iterate through to find rho / momentum
	for (int v = 0; v < nVels; v++) {
		rho[id] += f[id * nVels + v];
		u[id * dims + eX] += c[v][eX] * f[id * nVels + v];
		u[id * dims + eY] += c[v][eY] * f[id * nVels + v];
	}

	// Divide momentum by rho to calculate velocity
	u[id * dims + eX] = u[id * dims + eX] / rho[id];
	u[id * dims + eY] = u[id * dims + eY] / rho[id];
}

void GridClass::applyBC(int i, int j, int id) {

	eDirectionType normalDirection;
	std::vector<int> normalVector = getNormalVector(i, j, normalDirection);

	switch(type[id]) {

		case eWall:

			u_n[id * dims + eX] = 0.0;
			u_n[id * dims + eY] = 0.0;

			regularisedBC(i, j, id, normalVector, normalDirection);
			break;

		case eConvective:

			convectiveBC(j, id);
			break;

		case eFluid:

			break;

		case eVelocity:
			// Not implemented yet...
			break;
	}
}

void GridClass::regularisedBC(int i, int j, int id, std::vector<int> &normalVector, eDirectionType normalDirection) {
}

std::vector<int> GridClass::getNormalVector(int i, int j, eDirectionType &normalDirection) {

	std::vector<int> normalVector(dims, 0);

	// Get normal vector in x
	if (i == 0) {
		normalVector[eX] = 1;
		normalDirection = eX;
	}
	else if (i == Nx - 1) {
		normalVector[eX] = -1;
		normalDirection = eX;
	}

	// Get normal vector in y
	if (j == 0) {
		normalVector[eY] = 1;
		normalDirection = eY;
	}
	else if (j == Ny - 1) {
		normalVector[eY] = -1;
		normalDirection = eY;
	}

	// Check for corners
	if (normalVector[eX] != 0 && normalVector[eY] != 0) {

		if (type[(i + normalVector[eX]) * Ny + j] == eFluid) {
			normalVector[eY] = 0; 
			normalDirection = eX;
		}
		else if (type[i * Ny + j + normalVector[eY]] == eFluid) {
			normalVector[eX] = 0;
			normalDirection = eY;
		}
	}

	return normalVector;
}
