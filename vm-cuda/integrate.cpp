#include <iostream>
#include "vm.hpp"
#include "options.hpp"
#include <boost/foreach.hpp>
#include <omp.h>

int main(int argc, char *argv[]) {
    lamb_oseen_options ops(argc, argv);

    std::vector<particle> particles;
    read_particles(particles);

    float core_size = 2 * ops.h;
    VortexMethod vm(particles, core_size);

    float t = ops.t0, time_step = 0.01;
    for (unsigned iteration = 0; iteration < 100; ++iteration) {
        std::cout << "# iteration " << iteration << std::endl;
        double start = omp_get_wtime();
        vm.evaluate_velocity();
        vm.convect(time_step);
        t += time_step;
        double time = omp_get_wtime() - start;
        std::cout << "#   time measured: " << time << std::endl;
    }

}

