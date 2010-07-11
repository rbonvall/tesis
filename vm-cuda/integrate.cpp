#include <iostream>
#include "vm.hpp"
#include "options.hpp"
#include <boost/foreach.hpp>
#include <omp.h>
#if defined(CUDA_INTEGRATION)
#    include <cuda.h>
#endif

int main(int argc, char *argv[]) {
    lamb_oseen_options ops(argc, argv);

    std::vector<particle> particles;
    read_particles(particles);

    float core_size = 2 * ops.h;
    VortexMethod vm(particles, core_size);

    float t = ops.t0, time_step = 0.01;
    for (unsigned iteration = 0; iteration < 100; ++iteration) {
        std::cout << "##########################" << std::endl;
        std::cout << "### ITERATION " << iteration << std::endl;
        std::cout << "# t =  " << t << std::endl;

        double start = omp_get_wtime();
        vm.evaluate_velocity();
        vm.convect(time_step);
        double time = omp_get_wtime() - start;

        std::cout << "# iteration timing: " << time << std::endl;
        std::copy(particles.begin(), particles.end(),
                  std::ostream_iterator<particle>(std::cout, "\n"));
        std::cout << std::endl;

        t += time_step;
    }

}

