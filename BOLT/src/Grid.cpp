// Includes
#include "../include/Grid.h"
#include "../include/Utils.h"

// Main solver
void GridClass::solver() {

    // Run Grid Kernel
    lbmKernel();

    // Run Object Kernel
    if (objectPtr->hasIBM == true) {
        objectPtr->objectKernel();
    }
}

// Main LBM Kernel
void GridClass::lbmKernel() {
    
    if (WALL_RIGHT == eConvective) {
        convectiveSpeed();
    }

    // Set values at start of tStep
    f_n.swap(f);
    u_n.swap(u);
    rho_n.swap(rho);

    // Iterate through all points
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {

            int id = i * Ny + j;

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

// Streaming/colliding in one step - Pull Algorithm
void GridClass::streamCollide(int i, int j, int id) {

    for (int v = 0; v < nVels; v++) {
        int src_id = ((i - c[v][eX] + Nx) % Nx) * Ny + ((j - c[v][eY] + Ny) % Ny);
        f[id * nVels + v] = f_n[src_id * nVels + v] + omegaLocal[src_id] * (equilibrium(src_id, v) - f_n[src_id * nVels + v]) + latticeForce(src_id, v);
    }
}

// Equilibrium distribution function
double GridClass::equilibrium(int id, int v) {

    double C1 = SQ(u_n[id * dims + eX]) + SQ(u_n[id * dims + eY]);
    double C2 = u_n[id * dims + eX] * c[v][eX] + u_n[id * dims + eY] * c[v][eY];

    return rho_n[id] * w[v] * (1.0 + 3.0 * C2 + 4.5 * SQ(C2) - 1.5 * C1);
}

// Compute macroscopic properties
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
    u[id * dims + eX] = (u[id * dims + eX] + 0.5 * force_xyz[id * dims + eX]) / rho[id];
    u[id * dims + eY] = (u[id * dims + eY] + 0.5 * force_xyz[id * dims + eY]) / rho[id];
}

// Calculate lattice force
double GridClass::latticeForce(int id, int v) {

    return (1.0 - 0.5 * omegaLocal[id]) * (w[v] / SQ(c_s)) 
            * ((force_xyz[id * dims + eX] + force_ibm[id * dims + eX]) * (c[v][eX] - u_n[id * dims + eX] + c[v][eX] * (c[v][eX] * u_n[id * dims + eX] + c[v][eY] * u_n[id * dims + eY]) / SQ(c_s))
            + ((force_xyz[id * dims + eY] + force_ibm[id * dims + eY]) * (c[v][eY] - u_n[id * dims + eY] + c[v][eY] * (c[v][eX] * u_n[id * dims + eX] + c[v][eY] * u_n[id * dims + eY]) / SQ(c_s))));
}

// Apply boundary conditions
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

// Calculate convective boundary conditions
void GridClass::convectiveBC(int j, int id) {

    // Setting the f values assuming right wall is eConvective
    f[id * nVels + 2] = f_n[id * nVels + 2] + 3.0 * w[2] * (delU[j][eX] * c[2][eX] + delU[j][eY] * c[2][eY]);
    f[id * nVels + 6] = f_n[id * nVels + 6] + 3.0 * w[6] * (delU[j][eX] * c[6][eX] + delU[j][eY] * c[6][eY]);
    f[id * nVels + 8] = f_n[id * nVels + 8] + 3.0 * w[8] * (delU[j][eX] * c[8][eX] + delU[j][eY] * c[8][eY]);
}

// Calculate convective speed
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

// Regularise boundary conditions
void GridClass::regularisedBC(int i, int j, int id, std::vector<int> &normalVector, eDirectionType normalDirection) {

    // Implementation of method to regularise the boundary conditions - taking corners into account

    // Extrapolate if corner
    if (normalVector[eX] != 0 && normalVector[eY] != 0) {

        if (type[id] == eVelocity || type[id] == eWall)
            rho_n[id] = Utils::extrapolate(rho, normalVector, 1, i, j);
    }

    // If normal edge
    else {

        double fplus = 0.0, fzero = 0.0;

        // Loop through velocities
        for (int v = 0; v < nVels; v++) {

            // If normal is opposite then add to fplus
            if (c[v][normalDirection] == -normalVector[normalDirection]) {
                fplus += f[id * nVels + v];
            }

            // If it is perpendicular to wall then add to fzero
            else if (c[v][normalDirection] == 0) {
                fzero += f[id * nVels + v];
            }
        }

        // Velocity condition
        if (type[id] == eVelocity || type[id] == eWall) {
            rho_n[id] = (2.0 * fplus + fzero) / (1.0 - normalVector[normalDirection] * u_n[id * dims + normalDirection]);
        }
    }

    // Declare stresses
    double Sxx = 0.0, Syy = 0.0, Sxy = 0.0;

    // Update f values
    for (int v = 0; v < nVels; v++) {

        // Get feq
        double feq = equilibrium(id, v);

        // If corner then unknowns share at least one of normal components
        if (normalVector[eX] != 0 && normalVector[eY] != 0) {
            if (c[v][eX] == normalVector[eX] || c[v][eY] == normalVector[eY]) {

                // If buried link just set to feq
                if (std::inner_product(normalVector.begin(), normalVector.end(), c[v].begin(), 0) == 0) {
                    f[id * nVels + v] = feq;
                }
                else {
                    f[id * nVels + v] = feq + (f[id * nVels + opposite[v]] - equilibrium(id, opposite[v]));
                }
            }
        }

        // If other then unknowns share the normal vector component
        else {
            if (c[v][normalDirection] == normalVector[normalDirection]) {
                f[id * nVels + v] = feq + (f[id * nVels + opposite[v]] - equilibrium(id, opposite[v]));
            }
        }

        // Store off-equilibrium
        double fneq = f[id * nVels + v] - feq;

        // Compute off-equilbrium stress components
        Sxx += c[v][eX] * c[v][eX] * fneq;
        Syy += c[v][eY] * c[v][eY] * fneq;
        Sxy += c[v][eX] * c[v][eY] * fneq;
    }

    // Compute regularised non-equilibrium components and add to feq
    for (int v = 0; v < nVels; v++) {
        f[id * nVels + v] = equilibrium(id, v) + (w[v] / (2.0 * QU(c_s))) * (((SQ(c[v][eX]) - SQ(c_s)) * Sxx) + ((SQ(c[v][eY]) - SQ(c_s)) * Syy) + (2.0 * c[v][eX] * c[v][eY] * Sxy));
    }
}

// Get normal vector for specified lattice site
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

// Write information to terminal
void GridClass::writeInfo() {

    // Initialise maxVel to zero
    double maxVel = 0.0;

    // Iterate over all points to determine maxVel
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {

            int id = i * Ny + j;

            // Calculate basic vector magnitude
            double vel = std::sqrt(SQ(u[id * dims + eX]) + SQ(u[id * dims + eY]));

            if (vel > maxVel) {
                maxVel = vel;
            }
        }
    }

    // Calculate maxRe
    double maxRe = (maxVel * Dx / Dt) * ref_L / ref_nu;

    // Calculate timings
    loopTime = std::chrono::steady_clock::now() - startTime;
    loopTime /= t;

    // Output values
    std::cout << std::endl << std::endl;
    
    std::cout << "Time step " << t << " of " << nSteps << std::endl;
    std::cout << "Simulation has done " << t * Dt << " of " << nSteps * Dt << " seconds" << std::endl;
    std::cout << "MLUPS = " << Nx * Ny / (1000000.0 * loopTime.count()) << std::endl;

    std::cout << "Max Velocity = " << maxVel << std::endl;
    std::cout << "Max Velocity (m/s) = " << maxVel * Dx / Dt << std::endl;
    std::cout << "Max Reynolds Number = " << maxRe << std::endl;
}

// Write information to VTK file
void GridClass::writeVTK() {

    // Create file
    std::ofstream output;
    output.precision(PRECISION);
    output.open("Results/VTK/Fluid." + std::to_string(t) + ".vtk");

    // Write VTK header
    output << "# vtk DataFile Version 3.0\n";
    output << "FLUID\n";
    output << "ASCII\n";
    output << "DATASET RECTILINEAR_GRID\n";
    output << "DIMENSIONS " << Nx << " " << Ny << " 1\n";

    // Write out x-coordinates
    output << "X_COORDINATES " << Nx << " double\n";
    for (int i = 0; i < Nx; i++)
        output << i * Dx << " ";

    output << std::endl;

    // Write out y-coordinates
    output << "Y_COORDINATES " << Ny << " double\n";
    for (int j = 0; j < Ny; j++)
        output << j * Dx << " ";

    output << std::endl;

    // Write out z-coordinates
    output << "Z_COORDINATES " << 1 << " double\n";
    output << 0;

    output << std::endl;

    // Point data
    output << "POINT_DATA " << Nx * Ny << "\n";

    // Write density
    output << "SCALARS Density double 1\n";
    output << "LOOKUP_TABLE default\n";
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            output << rho[i * Ny + j] * Drho << "\n";
        }
    }

    // Write pressure
    output << "SCALARS Pressure double 1\n";
    output << "LOOKUP_TABLE default\n";
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            output << ref_P + (rho[i * Ny + j] - rho_p / Drho) * SQ(c_s) * Dm / (Dx * SQ(Dt)) << "\n";
        }
    }

    // Write velocity
    output << "VECTORS Velocity double\n";
    for (int j = 0; j < Ny; j++) {
        for (int i = 0; i < Nx; i++) {
            output << u[(i * Ny + j) * dims + eX] * (Dx / Dt) << " " << u[(i * Ny + j) * dims + eY] * (Dx / Dt) << " 0\n";
        }
    }

    // Close file
    output.close();
}

// Initialise grid values
void GridClass::initialiseGrid() {

#ifdef BLOCK
    int blockIdx = round((Ny - 1) * BLOCK) + 1;
    double blockH = height_p * BLOCK;
#else
    int blockIdx = 0;
    double blockH = 0.0;
#endif

    // Iterate through all points - set types
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {

            int id = i * Ny + j;

            // Left wall
            if (i == 0) {

                if (j < blockIdx) {
                    type[id] = eWall;   
                }
                else {
                    type[id] = WALL_LEFT;
                }
                
            }

            // Right wall
            if (i == Nx - 1) {

                if (j < blockIdx) {
                    type[id] = eWall;
                }
                else {
                    type[id] = WALL_RIGHT;
                }
            }

            // Bottom wall
            if (j == 0) {
                type[id] = WALL_BOTTOM;
            }

            // Top wall
            if (j == Ny - 1) {
                type[id] = WALL_TOP;
            }

            if (type[id] != eFluid) {
                BCVec.push_back(id);
            }
        }
    }

    // Set inlet velocity
    for (int j = 0; j < Ny; j++) {

#ifndef PROFILE

        u_in[j * dims + eX] = ux0_p * Dt / Dx;
        u_in[j * dims + eY] = uy0_p * Dt / Dx;

#else
        if (PROFILE == eParabolic) {

            double R = (height_p - blockH) / 2.0;
            double yPos = j * Dx - blockH - R;

            u_in[j * dims + eX] = 1.5 * (ux0_p * Dt / Dx) * (1.0 - SQ(yPos / R));
            u_in[j * dims + eY] = 1.5 * (uy0_p * Dt / Dx) * (1.0 - SQ(yPos / R));
        }
        else if (PROFILE == eBoundaryLayer) {
            std::cout << "eBoundaryLayer not yet supported..." << std::endl;
            exit(-1);
        }
#endif
    }


    // Iterate through all points - set initial velocity
    for (int i = 0; i < Nx; i++) {
        for (int j = 0; j < Ny; j++) {

            int id = i * Ny + j;

#ifdef PROFILE

            u[id * dims + eX] = u_in[j * dims + eX];
            u[id * dims + eY] = u_in[j * dims + eY];

#else

            u[id * dims + eX] = ux0_p * Dt / Dx;
            u[id * dims + eY] = uy0_p * Dt / Dx;

#endif

            if (type[id] == eWall) {
                u[id * dims + eX] = 0.0;
                u[id * dims + eY] = 0.0;
            }
        }
    }

    // Set values at beginning of tStep
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

// Default constructor
GridClass::GridClass() {

    // Output message
    std::cout << "Initialising GridClass." << std::endl;

    objectPtr = NULL;

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

    // Set timings
    startTime = std::chrono::steady_clock::now();
    // loopTime = 0.0;

    // Set array size and initialise
    u.resize(Nx * Ny * dims, 0.0);
    u_n.resize(Nx * Ny * dims, 0.0);
    rho.resize(Nx * Ny, rho0);
    rho_n.resize(Nx * Ny, rho0);

    omegaLocal.resize(Nx * Ny, omega);
    force_xyz.resize(Nx * Ny * dims, 0.0);
    force_ibm.resize(Nx * Ny * dims, 0.0);

    f.resize(Nx * Ny * nVels, 0.0);
    f_n.resize(Nx * Ny * nVels, 0.0);
    type.resize(Nx * Ny, eFluid);

    // Assuming right wall is eConvective for now...
    delU.resize(Ny, std::vector<double>(dims, 0.0));

    u_in.resize(Ny * dims, 0.0);

    initialiseGrid();
}
