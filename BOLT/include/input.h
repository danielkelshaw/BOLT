#ifndef INPUT_H
#define INPUT_H

// Include relevant headers
#include <omp.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <numeric>
#include <chrono>
#include <boost/filesystem.hpp>

#define THREADS 4
#define ORDERED

// Set resolution
const int resolution = 1;

// Lattice domain
const int Nx = 256;     // Number of lattice sites in x-direction
const int Ny = 256;     // Number of lattice sites in y-direction

// Physical domain
const double height_p = 1.0;        // Height of domain (m)
const double rho_p = 1.0;           // Fluid density (kg/m^3)
const double nu_p = 1.5111e-2;      // Fluid viscosity (m^2/s)

#define WALL_LEFT    eVelocity      // Boundary condition at left wall
#define WALL_RIGHT   eConvective    // Boundary condition at right wall
#define WALL_BOTTOM  eVelocity      // Boundary condition at bottom wall
#define WALL_TOP     eVelocity      // Boundary condition at top wall

// #define PROFILE eParabolic          // Inlet velocity profile
// #define BLOCK (1.0 / 3.0)           // Blockage ratio

// Initial conditions
const double ux0_p = 5.0;           // Initial x-velocity (m/s)
const double uy0_p = 0.0;           // Initial y-veloicty (m/s)

// Set time step
const double tStep = 5e-5 / (resolution * resolution);
const double omega = 1.0 / (nu_p * tStep / (std::pow(1.0 / std::sqrt(3.0), 2.0) * std::pow(height_p / (Ny - 1), 2.0)) + 0.5);

// Set simulation time
const double tSim = 4.0;
const int nSteps = static_cast<int>(std::round(tSim / tStep));

// Set output frequency
const int tinfo = nSteps / 10000;
const int tVTK = nSteps / 1000;

// Reference values
const double ref_nu = nu_p;
const double ref_rho = rho_p;
const double ref_P = 100000.0;
const double ref_L = 0.057;
const double ref_U = ux0_p;

// Set dimensions
const int dims = 2;

// Precision for outputs
const int PRECISION = 10;

// Enumerations
enum eDirectionType {eX, eY};
enum eLatType {eFluid, eWall, eVelocity, eConvective};
enum eProfile {eParabolic, eBoundaryLayer};
enum eBodyType {eCircle};

// Macros
#define SQ(x) ((x) * (x))
#define CU(x) ((x) * (x) * (x))
#define QU(x) ((x) * (x) * (x) * (x))

#endif