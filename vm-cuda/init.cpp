#include "vm.hpp"
#include "lamboseen.hpp"
#include "options.hpp"
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <boost/foreach.hpp>

#define EXIT_IF(cond, error_msg) \
        do if(cond) {std::cerr << error_msg << std::endl; return 1;} while (0);

int main(int argc, char *argv[]) {
    lamb_oseen_options ops(argc, argv);

    float x0 = ops.x0, x1 = ops.x1;
    float y0 = ops.y0, y1 = ops.y1;
    float h = ops.h;

    unsigned nr_cells = unsigned((x1 - x0) / h) *
                        unsigned((y1 - y0) / h);

    EXIT_IF(x0 >= x1 || y0 >= y1, "Degenerate geometry.");
    EXIT_IF(nr_cells == 0, "No room for particles.");

    lamb_oseen_vortex v(ops.gamma0, ops.nu);
    std::vector<particle> particles;
    particles.reserve(nr_cells);

    for (float x = x0 + h/2; x <= x1 - h/2; x += h)
        for (float y = y0 + h/2; y <= y1 - h/2; y += h) {
            float vort = v(x, y, ops.t0);
            float circ = vort * h * h;
            if (circ > ops.circulation_threshold) {
                particle p(x, y, vort * h * h);
                particles.push_back(p);
            }
        }

    if (nr_cells != particles.size())
        std::cerr << "nr_cells         = " << nr_cells << ", " <<
                     "particles.size() = " << particles.size() <<
                     std::endl;

    std::cout << "# Lamb-Oseen vortex particle discretization" << std::endl;
    std::cout << "# x0 = " << x0 << std::endl;
    std::cout << "# x1 = " << x1 << std::endl;
    std::cout << "# y0 = " << y0 << std::endl;
    std::cout << "# y1 = " << y1 << std::endl;
    std::cout << "# t0 = " << ops.t0 << std::endl;
    std::cout << "# total-circulation = " << ops.gamma0 << std::endl;
    std::cout << "# viscosity = " << ops.nu << std::endl;
    std::cout << "# circulation-threshold = " << ops.circulation_threshold << std::endl;
    std::cout << "# cell-size = " << h << std::endl;
    std::cout << "# nr-particles = " << particles.size() << std::endl;
    std::cout << std::endl;
    std::copy(particles.begin(), particles.end(),
              std::ostream_iterator<particle>(std::cout, "\n"));
}

