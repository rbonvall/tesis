#include <iostream>
#include <boost/foreach.hpp>
#include "vm.hpp"
#include "lamboseen.hpp"
#include "options.hpp"

int main(int argc, char *argv[]) {
    lamb_oseen_options ops(argc, argv);

    std::vector<particle> particles;
    read_particles(particles);

    std::cout << "# Analytical Lamb-Oseen velocity evaluation" << std::endl;
    std::cout << "# total-circulation = " << ops.gamma0 << std::endl;
    std::cout << "# viscosity = " << ops.nu << std::endl;
    std::cout << "# t0 = " << ops.t0 << std::endl;
    lamb_oseen_vortex v(ops.gamma0, ops.nu);
    BOOST_FOREACH(particle& p, particles) {
        std::pair<float, float> velocity = v.velocity(p.x, p.y, ops.t0);
        std::cout << p.x << ' ' << p.y << ' ' << p.circ << ' ' <<
                     velocity.first  << ' ' <<
                     velocity.second << std::endl;
    }
}

