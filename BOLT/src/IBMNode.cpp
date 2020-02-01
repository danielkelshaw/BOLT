// Includes
#include "../include/IBMNode.h"

void findSupport() {

}

void computeDs() {

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
