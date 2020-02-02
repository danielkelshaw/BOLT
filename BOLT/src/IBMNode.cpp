// Includes
#include "../include/IBMNode.h"

void findSupport() {

}

void computeDs() {

	// Set current_Ds as an arbitrarily large value
	double current_Ds = 100.0;

	// Loop through each node
	for (size_t n = 0; n < ibmPointer->nodes.size(); n++) {
		
		// Don't check self
		if (ID != ibmPointer->nodes[n].ID) {

			double mag = (pos - ibmPointer->nodes[n].pos) / ibmPointer->objectPtr->gridPtr->Dx;

			if (mag < current_Ds) {
				current_Ds = mag;
			}
		}		
	}

	// set ds value to current_Ds
	ds = current_Ds;
}

void interpolate() {

}

void forceCalc() {

}

void spread() {

}

void updateMacroscopic() {
	
}

IBMNodeClass::IBMNodeClass(IBMBodyClass *ibmPointer, int nodeID, const std::vector<double> &position) {

	ibmPtr = ibmPointer;

	ID = nodeID;

	pos	= position;
	pos0 = position;

	vel.resize(dims, 0.00;)
	force.resize(dims, 0.00;)

	ds = 0.0;
	epsilon = 0.0;

	interp_mom.resize(dims, 0.0);
	interp_rho = 0.0;
}
