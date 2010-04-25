#include <iostream>
#include "vm.hpp"
#include <boost/foreach.hpp>

int main(int argc, char *argv[]) {
    std::vector<particle> particles;
    read_particles(particles);

    VortexMethod vm(particles, 3e-2);
    vm.evaluate_velocity();

    std::cout << "# Particles after velocity evaluation" << std::endl;
    BOOST_FOREACH(particle& p, vm.particles) {
        std::cout << p.x << ' ' << p.y << ' ' << p.circ << ' ' <<
                     p.u << ' ' << p.v << std::endl;
    }


}

