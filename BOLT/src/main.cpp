#include "../include/Grid.h"

int main() {
	std::cout << "BOLT: LBM Simulator" << std::endl;

	GridClass grid;

	for (grid.t = 0; grid.t <= nSteps; grid.t++) {
		grid.solver();

		if (grid.t % 250 == 0) {
			std::cout << "Step: " << grid.t << "/" << nSteps << std::endl;
		}
	}
}
