#ifndef IBMBODY_H
#define IBMBODY_H

// Includes
#include "input.h"
#include "Objects.h"

// forward declarations
class ObjectsClass;

// IBMBody Class
class IBMBodyClass {

public:

	// defaul constructor / destructor
	IBMBodyClass() {};
	~IBMBodyClass() {};

	// custom constructor for circle case
	IBMBodyClass(ObjectsClass *objects, int bodyID, const std::vector<double> &pos, double radius);

public:

private:

	ObjectsClass *objectPtr;

	int ID;
	eBodyType bodyType;
};


#endif