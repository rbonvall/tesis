#pragma once

#include <vector>
#include <iostream>

struct particle {
    float x, y; // position
    float circ; // circulation
    float u, v; // velocity

    particle(float x, float y, float circ) :
        x(x), y(y), circ(circ), u(0), v(0) {}

};

std::ostream& operator <<(std::ostream& out, particle p)  {
    out << p.x << ' ' << p.y << ' ' <<
           p.circ << ' ' << p.u << ' ' << p.v;
    return out;
}


class VortexMethod {
    public:
    void start(unsigned nr_iterations);

    VortexMethod(std::vector<particle>& particles) :
        particles(particles) {}

    private:
    std::vector<particle> particles;

    void evaluate_velocity();
    void convect();
};



void read_particles(std::vector<particle>& particles, std::istream& in = std::cin) {
    while (in) {
        float x, y, circ;
        in >> x >> y >> circ;
        particles.push_back(particle(x, y, circ));
    }
}
