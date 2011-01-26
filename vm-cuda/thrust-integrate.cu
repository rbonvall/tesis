#include <vector>
#include <algorithm>
#include "vm.hpp"
#include "thrust/host_vector.h"
#include "thrust/device_vector.h"
#include "thrust/transform.h"


struct particle_to_float4 {
    float4 operator()(particle p) {
        return make_float4(p.x, p.y, p.circ, 0.0);
    }
};


__global__ void
integrate(float dt, unsigned nr_particles, float4 *old_part, float4 *new_part) {
    unsigned pid = blockIdx.x * blockDim.x + threadIdx.x;

    // fetch particle from global memory
    float4 p = old_part[pid];

    // compute velocity and derivative of circulation for particle p
    float3 derivatives /*= eval_derivatives(p, old_part, nr_particles)*/;

    // compute the derivative of the circulation by using PSE

    // convect particle and copy it to global memory
    p.x += derivatives.x * dt;
    p.y += derivatives.y * dt;
    p.z += derivatives.z * dt;
    new_part[pid] = p;
}


void solve(std::vector<particle> particles, float dt, unsigned nr_iterations) {
    unsigned N = particles.size();

    thrust::host_vector<float4> ps_h(N);
    thrust::transform(particles.begin(), particles.end(),
                      ps_h.begin(),
                      particle_to_float4());

    thrust::device_vector<float4> ps_d[2] = { // I hate you, C++.
        thrust::device_vector<float4>(N),
        thrust::device_vector<float4>(N),
    };
    unsigned current_read = 0, current_write = 1;

    thrust::copy(ps_h.begin(), ps_h.end(), ps_d[current_read].begin());

    for (unsigned i = 0; i < nr_iterations; ++i) {
        integrate<<<N/256, 256>>>(dt, N,
                (float4*) thrust::raw_pointer_cast(&ps_d[current_read]),
                (float4*) thrust::raw_pointer_cast(&ps_d[current_write]));
        cudaSynchronizeThreads();

        std::swap(current_read, current_write);
    }

}

