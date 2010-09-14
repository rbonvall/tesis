#pragma once

#include <vector>
#include <iostream>
#include <functional>
#include <cmath>
#include <limits>

struct particle {
    float x, y; // position
    float circ; // circulation
    float u, v; // velocity

    particle(float x, float y, float circ) :
        x(x), y(y), circ(circ), u(0), v(0) {}

};

std::ostream& operator <<(std::ostream& out, particle p);
void read_particles(std::vector<particle>& particles, std::istream& in = std::cin);

//struct squared_distance :
//    public std::binary_function<const particle&, const particle&, float>
//{
//    float operator() (const particle& p, const particle& q) {
//        float dx = p.x - q.x;
//        float dy = p.y - q.y;
//        return dx * dx + dy * dy;
//    }
//};

struct squared_distance :
    public std::binary_function<const particle&, const particle&, float>
{
    float operator() (const particle& p, const particle& q) {
        float dx = p.x - q.x;
        float dy = p.y - q.y;
        return dx * dx + dy * dy;
    }
};


class VortexMethod {
    public:
    void start(unsigned nr_iterations);

    VortexMethod(std::vector<particle>& particles, float core_size) :
        particles(particles), core_size(core_size) {}

    std::vector<particle> particles;
    const float core_size;

    void evaluate_velocity();
    void convect(float time_step);
    void diffuse(float viscosity);
};




struct biot_savart_kernel_factor :
    public std::unary_function<float, float>
{};

struct gaussian_bskf :
    public biot_savart_kernel_factor
{
    float sq_epsilon;

    gaussian_bskf(float epsilon) : sq_epsilon(epsilon * epsilon) {}

    float operator() (float sq_r) {
        return (1 - exp(-sq_r/sq_epsilon)) /
               (2 * M_PI * (sq_r + std::numeric_limits<float>::epsilon()));
    }
};

