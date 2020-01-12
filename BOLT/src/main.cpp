#include "../include/Grid.h"

int main() {
	std::cout << "BOLT: LBM Simulator" << std::endl;

	GridClass grid;

	grid.solver();

	std::cout << "Iteration run succesfully..." << std::endl;
}
