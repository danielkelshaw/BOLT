#include "../include/Grid.h"
#include "../include/Objects.h"

int main() {
    std::cout << "BOLT: LBM Simulator" << std::endl;

    GridClass grid;
    ObjectsClass objects(grid);

    for (grid.t = 1; grid.t <= nSteps; grid.t++) {
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
