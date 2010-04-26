#include "vm.hpp"
#include "lamboseen.hpp"
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include <boost/foreach.hpp>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#define EXIT_IF(cond, error_msg) \
        do if(cond) {std::cerr << error_msg << std::endl; return 1;} while (0);

int main(int argc, char *argv[]) {
    float gamma0, nu, x0, x1, y0, y1, h, t0, circulation_threshold;

    po::options_description desc;
#   define OPTION(op, var, defval) ((op), po::value<float>(&(var))->default_value(defval))
    desc.add_options()
        OPTION("total-circulation", gamma0, 1.0)
        OPTION("viscosity",         nu,     5e-4)
        OPTION("circ-threshold",    circulation_threshold, 1e-5)
        OPTION("cell-size",         h,      7.8125e-3)
        OPTION("x0", x0, -0.3)
        OPTION("x1", x1, +0.3)
        OPTION("y0", y0, -0.3)
        OPTION("y1", y1, +0.3)
        OPTION("t",  t0, 4.00)
    ;
#   undef OPTION
    po::variables_map vars;
    po::store(po::parse_command_line(argc, argv, desc), vars);
    po::notify(vars);

    unsigned nr_cells = static_cast<unsigned>((x1 - x0) / h) *
                        static_cast<unsigned>((y1 - y0) / h);

    EXIT_IF(x0 >= x1 || y0 >= y1, "Degenerate geometry.");
    EXIT_IF(nr_cells == 0, "No room for particles.");

    lamb_oseen_vortex v(gamma0, nu);
    std::vector<particle> particles;
    particles.reserve(nr_cells);

    for (float x = x0 + h/2; x <= x1 - h/2; x += h)
        for (float y = y0 + h/2; y <= y1 - h/2; y += h) {
            float vort = v(x, y, t0);
            float circ = vort * h * h;
            if (circ > circulation_threshold) {
                particle p(x, y, vort * h * h);
                particles.push_back(p);
            }
        }

    if (nr_cells != particles.size())
        std::cerr << "nr_cells         = " << nr_cells << ", " <<
                     "particles.size() = " << particles.size() <<
                     std::endl;

    std::copy(particles.begin(), particles.end(),
              std::ostream_iterator<particle>(std::cout, "\n"));
}

