#include <vector>
#include <algorithm>
#include "vm.hpp"
#include "thrust/host_vector.h"
#include "thrust/device_vector.h"
#include "thrust/transform.h"
#include <cmath>

#define T 512 // threads per block

#define SQ(x) ((x) * (x))
const float ONE_OVER_TWO_PI = (1 / (2 * M_PI));
const float E = 9e-2;
const float SQ_E = SQ(E);
const float ONE_OVER_SQ_E = 1 / SQ_E;
const float VISCOSITY = 5e-3;

struct particle_to_float4 {
    float4 operator()(particle p) {
        return make_float4(p.x, p.y, p.circ, 0.0);
    }
};


__device__ float eta(float sq_r) {
    return ONE_OVER_TWO_PI * expf(-sq_r/SQ_E);
}

__device__ float k_factor(float sq_r) {
    return (ONE_OVER_TWO_PI / sq_r) * (1 - expf(-sq_r / SQ_E));
}


__device__ float3
particle_interaction(float3 d, float4 p, float4 q) {

    float2 r = make_float2(p.x - q.x, p.y - q.y);
    float sq_r = SQ(r.x) + SQ(r.y);

    // update velocity
    float vel_kernel_factor = k_factor(sq_r);
    d.x += vel_kernel_factor * -r.y;
    d.y += vel_kernel_factor *  r.x;

    // update circulation
    d.z += (q.z - p.z) * eta(sq_r);

    return d;
}


__device__ float3
update_tile(float4 p, float3 d) {
    extern __shared__ float4 shared_particles[];

    unsigned long i = 0;
    unsigned int counter = 0;

#   define SHARED(i) (shared_particles[(i) + blockDim.x * threadIdx.x])

    while (counter < blockDim.x) {
        d = particle_interaction(d, p, SHARED(i++));
        ++counter;
    }

    d.z *= VISCOSITY * ONE_OVER_SQ_E;
    return d;
}


__device__ float3
eval_derivatives(float4 p, float4 *particles, float N) {
    unsigned pid = blockIdx.x * blockDim.x + threadIdx.x;

    unsigned num_tiles = 0; // TODO

    extern __shared__ float4 shared_particles[];

    float3 derivatives = make_float3(0.0f, 0.0f, 0.0f);
    for (int tile = 0; tile < num_tiles; ++tile) {
        shared_particles[threadIdx.x] = particles[0]; // TODO
        __syncthreads();

        derivatives = update_tile(p, derivatives);
        __syncthreads();
    }

    return derivatives;
}


__global__ void
integrate(float dt, unsigned nr_particles, float4 *old_particles, float4 *new_particles) {
    unsigned pid = blockIdx.x * blockDim.x + threadIdx.x;

    // fetch particle from global memory
    if (pid < nr_particles) {
        float4 p = old_particles[pid];

        // compute velocity and derivative of circulation for particle p
        float3 derivatives = eval_derivatives(p, old_particles, nr_particles);

        // convect particle and copy it to global memory
        p.x += derivatives.x * dt;
        p.y += derivatives.y * dt;
        p.z += derivatives.z * dt;
        new_particles[pid] = p;
    }
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
        integrate<<<std::ceil(N / T), T, T * sizeof(float4)>>>(dt, N,
                (float4*) thrust::raw_pointer_cast(&ps_d[current_read]),
                (float4*) thrust::raw_pointer_cast(&ps_d[current_write]));
        cudaThreadSynchronize();

        std::swap(current_read, current_write);
    }

}

