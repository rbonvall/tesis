#pragma once

#include <cstdlib>
#include <boost/program_options.hpp>
namespace po = boost::program_options;

struct lamb_oseen_options {
    float gamma0, nu, x0, x1, y0, y1, h, t0, circulation_threshold;

    lamb_oseen_options(int argc, char** argv) {
        po::options_description desc;
#       define OPTION(op, var, defval) ((op), po::value<float>(&(var))->default_value(defval))
        desc.add_options()
            ("help", "Print a help message")
            OPTION("total-circulation", gamma0, 1.0)
            OPTION("circ-threshold", circulation_threshold, 1e-5)
            OPTION("viscosity",   nu, 5e-4)
            OPTION("cell-size,h", h,  7.8125e-3)
            OPTION("x0", x0, -0.3)
            OPTION("x1", x1, +0.3)
            OPTION("y0", y0, -0.3)
            OPTION("y1", y1, +0.3)
            OPTION("t",  t0, 4.00)
        ;
#       undef OPTION
        po::variables_map vars;
        po::store(po::parse_command_line(argc, argv, desc), vars);
        po::notify(vars);

        if (vars.count("help")) {
            std::cout << desc << std::endl;
            std::exit(0);
        }
    }
};
