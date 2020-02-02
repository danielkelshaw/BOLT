#ifndef IBMBODY_H
#define IBMBODY_H

// Includes
#include "input.h"
#include "Objects.h"
#include "IBMNode.h"

// forward declarations
class ObjectsClass;
class IBMNodeClass

// IBMBody Class
class IBMBodyClass {

// constructors / destructors
public:

	// default constructor / destructor
	IBMBodyClass() {};
	~IBMBodyClass() {};

	// custom constructor for circle case
	IBMBodyClass(ObjectsClass *objects, int bodyID, const std::vector<double> &pos, double radius);

// public members
public:

// private members
private:

	ObjectsClass *objectPtr;

	int ID;
	eBodyType bodyType;

	std::vector<IBMNodeClass> nodes;
};


#endif
