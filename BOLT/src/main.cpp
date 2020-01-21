#include "../include/Grid.h"

int main() {
    std::cout << "BOLT: LBM Simulator" << std::endl;

    GridClass grid;

    for (grid.t = 0; grid.t <= nSteps; grid.t++) {
        grid.solver();

        // Write out information every tinfo steps
        if (grid.t % tinfo == 0) {
            grid.writeInfo();
        }

        // Write VTK every tVTK steps
        if (grid.t % tVTK == 0) {
            grid.writeVTK();
        }
    }
}
