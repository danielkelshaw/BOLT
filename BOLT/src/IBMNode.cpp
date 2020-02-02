// Includes
#include "../include/IBMNode.h"

void IBMNodeClass::findSupport() {

    // clear existing support points
    supps.clear();

    double stencilWidth = 1.5;

    // get lattice spacing
    double Dx = ibmPtr->objectPtr->gridPtr->Dx;

    // closest lattice sites to node
    int i_near = static_cast<int>(round(pos[eX] / Dx));
    int j_near = static_cast<int>(round(pos[eY] / Dx));

    for (int i = i_near - 2; i < i_near + 2; i++) {

        double distX = fabs(pos[eX] / Dx - i);

        for (int j = j_near - 2; j < j_near + 2; j++) {

            double distY = fabs(pos[eY] / Dx - j);

            if (distX < stencilWidth && distY < stencilWidth && i >= 0 && i <= Nx - 1 && j >= 0 && j <= Ny - 1) {
                supps.emplace_back(i, j, Utils::diracDelta(distX) * Utils::diracDelta(distY));
            }
        }
    }
}

void IBMNodeClass::computeDs() {

    // Set current_Ds as an arbitrarily large value
    double current_Ds = 10.0;

    // Loop through each node
    for (size_t n = 0; n < ibmPtr->nodes.size(); n++) {
        
        // Don't check self
        if (ID != ibmPtr->nodes[n].ID) {

            double mag = sqrt((pos - ibmPtr->nodes[n].pos) * (pos - ibmPtr->nodes[n].pos)) / ibmPtr->objectPtr->gridPtr->Dx;

            if (mag < current_Ds) {
                current_Ds = mag;
            }
        }       
    }

    // set ds value to current_Ds
    ds = current_Ds;
}

void IBMNodeClass::interpolate() {

    // Get pointer to grid
    GridClass *gridPtr = ibmPtr->objectPtr->gridPtr;

    interp_rho = 0.0;
    fill(interp_mom.begin(), interp_mom.end(), 0.0);

    // Loop through all points and interpolate
    for (size_t s = 0; s < supps.size(); s++) {

        // Find ID of point
        int id = supps[s].idx * Ny + supps[s].jdx;

        // Interpolate density
        interp_rho += gridPtr->rho[id] * supps[s].diracVal;

        // Interpolate momentum
        for (int d = 0; d < dims; d++) {
            interp_mom[d] += gridPtr->rho[id] * gridPtr->u[id * dims + d] * supps[s].diracVal;
        }
    }
}

void IBMNodeClass::forceCalc() {

    double velScale = ibmPtr->objectPtr->gridPtr->Dt / ibmPtr->objectPtr->gridPtr->Dx;
    force = 2.0 * (velScale * interp_rho * vel - interp_mom);
}

void IBMNodeClass::spread() {

    GridClass *gridPtr = ibmPtr->objectPtr->gridPtr;

    for (size_t s = 0; s < supps.size(); s++) {

        int id = supps[s].idx * Ny + supps[s].jdx;

        double Fx = force[eX] * epsilon * ds * supps[s].diracVal;
        double Fy = force[eY] * epsilon * ds * supps[s].diracVal;

        gridPtr->force_ibm[id * dims + eX] += Fx;
        gridPtr->force_ibm[id * dims + eY] += Fy;
    }
}

void IBMNodeClass::updateMacroscopic() {

    GridClass *gridPtr = ibmPtr->objectPtr->gridPtr;
    int nVels = gridPtr->nVels;

    for (size_t s = 0; s < supps.size(); s++) {

        double rhoTmp = 0.0;
        double uTmp = 0.0;
        double vTmp = 0.0;

        int id = supps[s].idx * Ny + supps[s].jdx;

        for (int v = 0; v < nVels; v++) {
            rhoTmp += gridPtr->f[id * nVels + v];
            uTmp   += gridPtr->c[v][eX] * gridPtr->f[id * nVels + v];
            uTmp   += gridPtr->c[v][eY] * gridPtr->f[id * nVels + v];
        }

        uTmp = (uTmp + 0.5 * (gridPtr->force_xyz[id * dims + eX] + gridPtr->force_ibm[id * dims + eX])) / rhoTmp;
        vTmp = (vTmp + 0.5 * (gridPtr->force_xyz[id * dims + eY] + gridPtr->force_ibm[id * dims + eY])) / rhoTmp;

        gridPtr->rho[id] = rhoTmp;
        gridPtr->u[id * dims + eX] = uTmp;
        gridPtr->u[id * dims + eY] = vTmp;
    }   
}

IBMNodeClass::IBMNodeClass(IBMBodyClass *ibmPointer, int nodeID, const std::vector<double> &position) {

    ibmPtr = ibmPointer;

    ID = nodeID;

    pos = position;
    pos0 = position;

    vel.resize(dims, 0.0);
    force.resize(dims, 0.0);

    ds = 0.0;
    epsilon = 0.0;

    interp_mom.resize(dims, 0.0);
    interp_rho = 0.0;
}
