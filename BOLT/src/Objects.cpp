// Includes
#include "../include/Objects.h"

void ObjectsClass::readGeometry() {

	// open config file
	std::ifstream file;
	file.open("input/geometry.config", std::ios::in);

	// handle failure to open file
	if (!file.is_open()) {
		std::cout << "Error opening geometry config file..." << std::endl;
		exit(-1);
	}

	// skip comments in config file
	std::streamoff fileOffset;
	std::string line;
	file.seekg(std::ios::beg);

	do {
		fileOffset = file.tellg();
		std::getline(file, line);
	} while (line[0] == '#' && !file.eof());

	// move cursor over comments
	file.seekg(fileOffset, std::ios::beg);

	std::string bodyCase;

	// read in geometries
	int bodyID = 0;
	while (file) {

		// get type of body
		file >> bodyCase;

		int number;
		file >> number;

		std::vector<double> start(dims);
		file >> start[eX];
		file >> start[eY];

		std::vector<double> space(dims);
		file >> space[eX];
		file >> space[eY];

		if (bodyCase == "CIRCLE") {

			double radius;
			file >> radius;

			for (int i = 0; i < number; i++) {

				std::vector<double> pos = {start[eX] + i * space[eX], start[eY] + i * space[eY]};
				ibmBody.emplace_back(this, bodyID, pos, radius);
				bodyID++;
			}

		}
		else {
			std::cout << "Only CIRCLE is supported at the moment..." << std::endl;
			exit(-1);
		}

		bodyCase = "NONE";
	}

	if (ibmBody.size() > 0) {
		hasIBM = true;
	}

	if (hasIBM == true) {
		std::cout << "found " << ibmBody.size() << " object(s)" << std::endl;
	}
	else if (hasIBM == false) {
		std::cout << "no objects found" << std::endl;
	}

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

	hasIBM = false;

	// Read in geometry
	readGeometry();

	// Initialise objects
	initialiseObjects();
}
