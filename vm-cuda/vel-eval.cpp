#include <iostream>
#include "vm.hpp"
#include "options.hpp"
#include <boost/foreach.hpp>

int main(int argc, char *argv[]) {
    lamb_oseen_options ops(argc, argv);

    std::vector<particle> particles;
    read_particles(particles);

    float core_size = 3e-2;
    VortexMethod vm(particles, core_size);
    vm.evaluate_velocity();

    std::cout << "# Particles after velocity evaluation" << std::endl;
    BOOST_FOREACH(particle& p, vm.particles) {
        std::cout << p.x << ' ' << p.y << ' ' << p.circ << ' ' <<
                     p.u << ' ' << p.v << std::endl;
    }


}

