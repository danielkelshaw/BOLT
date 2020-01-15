#include "../include/Grid.h"

int main() {
	std::cout << "BOLT: LBM Simulator" << std::endl;

	GridClass grid;

	for (grid.t = 0; grid.t <= nSteps; grid.t++) {
		grid.solver();

		// Write out information every 250 tSteps
		if (grid.t % 250 == 0) {
			grid.writeInfo();
			grid.writeVTK();
		}
	}
}
