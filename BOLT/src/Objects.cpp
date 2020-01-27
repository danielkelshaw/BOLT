// Includes
#include "../include/Objects.h"

void ObjectsClass::readGeometry() {

}

void ObjectsClass::initialiseObjects() {
	
}

ObjectsClass::ObjectsClass(GridClass &grid) {

	std::cout << "Initialising Objects..." << std::endl;

	// Set pointers
	gridPtr = &grid;
	gridPtr->objectPtr = this;

	// Initial sub-iteration values
	subIt = 0;
	subRes = 0.0;
	relax = 1.0;
	subNum = 0.0;
	subDen = 0.0;

	// Read in geometry
	readGeometry();

	// Initialise objects
	initialiseObjects();
}
