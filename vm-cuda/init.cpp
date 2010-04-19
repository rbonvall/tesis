#include "vm.hpp"
#include "lamboseen.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <boost/foreach.hpp>

int main(int argc, char *argv[]) {
    float gamma0 = 1.0;
    float nu = 5e-4;
    lamb_oseen_vortex v(gamma0, nu);


    float x0 = -0.3, x1 = 0.3;
    float y0 = -0.3, y1 = 0.3;
    float h = 7.8125e-3;
    unsigned int nr_cells = static_cast<int>((x1 - x0) / h) *
                            static_cast<int>((y1 - y0) / h);

    std::vector<particle> particles;
    particles.reserve(nr_cells);

    float t0 = 4.00;
    float circulation_threshold = 1e-5;
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

