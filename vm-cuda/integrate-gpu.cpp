#include <iostream>
#include "vm.hpp"
#include "options.hpp"
#include "cuda-integrate.hpp"
#include <boost/foreach.hpp>
#include <omp.h>
#include <vector>

int main(int argc, char *argv[]) {
    lamb_oseen_options ops(argc, argv);

    //std::vector<float> x, y, circ, u, v;
    //float xp, yp, circp, up, vp;
    //while (std::cin >> xp >> yp >> circp >> up >> vp) {
    //    x.push_back(xp);
    //    y.push_back(yp);
    //    circ.push_back(circp);
    //    // ignore up and vp
    //}

    std::vector<particle> particles;
    read_particles(particles);
    float x, y, circ, u, v;

    gpu_init(particles);
    std::cout << "Particles copied to GPU" << std::endl;

    float core_size = 2 * ops.h;
    float t = ops.t0, time_step = 0.01;
    for (unsigned iteration = 0; iteration < 100; ++iteration) {
        std::cout << "##########################" << std::endl;
        std::cout << "### ITERATION " << iteration << std::endl;
        std::cout << "# t =  " << t << std::endl;

        double start = omp_get_wtime();
        vm_integrate(time_step, 1, 2);
        double time = omp_get_wtime() - start;

        t += time_step;
    }

    gpu_finalize();


}

