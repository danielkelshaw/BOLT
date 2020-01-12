#ifndef INPUT_H
#define INPUT_H

// Include relevant headers
#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>

// Set resolution
const int resolution = 1;

// Lattice domain
const int Nx = 256;     // Number of lattice sites in x-direction
const int Ny = 128;     // Number of lattice sites in y-direction

// Physical domain
const double height_p = 1.0;
const double rho_p = 1.0;           // Fluid density (kg/m^3)
const double nu_p = 1.5111e-2;      // Fluid viscosity (m^2/s)

// Initial conditions
const double ux0_p = 5.0;           // Initial x-velocity (m/s)
const double uy0_p = 0.0;           // Initial y-veloicty (m/s)


// Set time step
const double tStep = 5e-5 / (resolution * resolution);
const double omega = 1.0 / (nu_p * tStep / (std::pow(1.0 / std::sqrt(3.0), 2.0) * std::pow(height_p / (Ny - 1), 2.0)) + 0.5);

const double tSim = 4.0;
const int nSteps = static_cast<int>(std::round(tSim / tStep));

// Set dimensions
const int dims = 2;

// Enumerations
enum eDirectionType {eX, eY};
enum eLatType {eFluid, eWall, eVelocity, eConvective};

// Macros
#define SQ(x) ((x) * (x))
#define CU(x) ((x) * (x) * (x))
#define QU(x) ((x) * (x) * (x) * (x))

#endif