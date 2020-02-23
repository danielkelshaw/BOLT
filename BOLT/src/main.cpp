#include "../include/Grid.h"
#include "../include/Objects.h"

int main() {
    std::cout << "BOLT: LBM Simulator" << std::endl;

    double start = omp_get_wtime();

#ifdef THREADS
    omp_set_num_threads(THREADS);
#endif

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
            
            objects.writeVTK();
            objects.writeForces();
        }
    }

    std::cout << "Simulation Finished." << std::endl;

    double end = omp_get_wtime();
    std::cout << "Simulation took " << end - start << " secs" << std::endl;
}
