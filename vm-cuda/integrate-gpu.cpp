#include <iostream>
#include "vm.hpp"
#include "options.hpp"
#include "thrust-integrate.h"
#include <boost/foreach.hpp>
#include <omp.h>
#include <vector>

int main(int argc, char *argv[]) {
    lamb_oseen_options ops(argc, argv);

    std::vector<particle> particles;
    std::vector<particle> velocities;
    read_particles(particles);

    double start = omp_get_wtime();
    float time_step = 0.01;
    unsigned nr_iterations=100;
    solve(particles, time_step, nr_iterations);
    double time = omp_get_wtime() - start;

    std::cout << "Time step: " << time_step << std::endl;
    std::cout << "Nr.iterations: " << nr_iterations << std::endl;
    std::cout << "Nr.particles: " << particles.size() << std::endl;
    std::cout << "Time: " << time << " seconds" << std::endl;

}

