#ifndef IBMNODE_H
#define IBMNODE_H

#include "input.h"
#include "IBMBody.h"

// Forward declarations
class IBMBodyClass;

class IBMNodeClass {

public:
	IBMNodeClass(IBMBodyClass *ibmPtr, int nodeID, const std::vector<double> &pos);
	~IBMNodeClass() {};

private:

	IBMBodyClass *ibmPtr;

	int ID;

	std::vector<double> pos;         // position of the node
	std::vector<double> pos0;        // initial position of the node
	std::vector<double> vel;         // velocity at the node
	std::vector<double> force;       // force

	// spacing and force multiplier
	double ds;                       // spacing between points
	double epsilon;                  // multiplier for IBM

	// interpolated values
	std::vector<double> interp_mom;  // interpolated momentum
	double interp_rho                // interpolated density

private:

	// IBM methods
	void findSupport();              // find support
	void computeDs();                // compute spacing between nodes
	void interpolate();              // interpolation
	void forceCalc();                // force calculation
	void spread();                   // spread force back
	void updateMacroscopic()         // update macroscopic values at supports
};

#endif
