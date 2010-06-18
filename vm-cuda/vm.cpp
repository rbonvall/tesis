#include "vm.hpp"
#include <cmath>
#include <algorithm>
#include <functional>
#include <numeric>
#include <boost/foreach.hpp>
//#include <boost/bind.hpp>
//#include <boost/lambda/lambda.hpp>

std::ostream& operator <<(std::ostream& out, particle p)  {
    out << p.x << ' ' << p.y << ' ' <<
           p.circ << ' ' << p.u << ' ' << p.v;
    return out;
}

void read_particles(std::vector<particle>& particles, std::istream& in) {
    float x, y, circ, u, v;
    while (in >> x >> y >> circ >> u >> v)
        particles.push_back(particle(x, y, circ));
}


void VortexMethod::evaluate_velocity() {
    gaussian_bskf K(core_size);

    unsigned N = particles.size();
    std::vector<float> bs_kernel_factors(N);

#if 0
    BOOST_FOREACH(particle& p, particles) {
        // Compute the kernel factor K(xp - xq)
        // against all other particles q
        std::transform(particles.begin(), particles.end(),
                       bs_kernel_factors.begin(),
                       boost::bind(K, boost::bind(squared_distance(), p, ::_1)));

        std::inner_product
#endif
    BOOST_FOREACH(particle& p, particles) {
        p.u = p.v = 0;
    }

    for (unsigned i = 0; i < N; ++i) {
        // fetch i-th particle and reset its velocity
        particle& p(particles[i]);

        // compute mutual velocity contributions
        // against each other particle q with a
        // higher index (j > i)
        for (unsigned j = i + 1; j < N; ++j) {
            particle& q(particles[j]);
            float dx = p.x - q.x;
            float dy = p.y - q.y;
            float kf = K(dx * dx + dy * dy);

            p.u += q.circ * -dy * kf;
            p.v += q.circ *  dx * kf;

            q.u -= p.circ * -dy * kf;
            q.v -= p.circ *  dx * kf;

        }
    }
}

void VortexMethod::convect(float time_step) {
    BOOST_FOREACH(particle& p, particles) {
        p.x += p.u * time_step;
        p.y += p.v * time_step;
    }
}
void VortexMethod::diffuse(float viscosity) {
}

