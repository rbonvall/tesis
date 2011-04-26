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
update_tile(float4 p, float3 derivatives) {
    extern __shared__ float4 shared_particles[];

#   define SHARED(i) (shared_particles[(i) + blockDim.x * threadIdx.x])

    unsigned long i = 0;
    unsigned int counter = 0;
    while (counter < blockDim.x) {
        derivatives = particle_interaction(derivatives, p, SHARED(i++));
        ++counter;
    }

    derivatives.z *= VISCOSITY * ONE_OVER_SQ_E;
    return derivatives;
}


__device__ float3
eval_derivatives(float4 p, float4 *particles, float N) {
    extern __shared__ float4 shared_particles[];
    unsigned const t = threadIdx.x;
    unsigned num_tiles = std::ceil(N / T);

    float3 derivatives = make_float3(0.0f, 0.0f, 0.0f);
    for (int tile = 0; tile < num_tiles; ++tile) {
        shared_particles[t] = particles[tile * T + t];
        __syncthreads();

        derivatives = update_tile(p, derivatives);
        __syncthreads();
    }

    return derivatives;
}


template <bool get_derivatives>
__global__ void
integrate(float dt, unsigned nr_particles, float4 *old_particles, float4 *new_particles, float4 *new_derivatives) {
    unsigned pid = blockIdx.x * blockDim.x + threadIdx.x;

    // fetch particle from global memory
    if (pid < nr_particles) {
        float4 p = old_particles[pid];

        // compute velocity and derivative of circulation for particle p
        float3 derivatives = eval_derivatives(p, old_particles, nr_particles);

        // integrate trajectories and circulation
        p.x += derivatives.x * dt;
        p.y += derivatives.y * dt;
        p.z += derivatives.z * dt;

        // put the particle back in global memory
        new_particles[pid] = p;

        if (get_derivatives) {
            new_derivatives[pid] = make_float4(derivatives.x, derivatives.y, derivatives.z, 0.0f);
        }
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
        unsigned nr_blocks = static_cast<unsigned>(std::ceil(N / T));
        integrate<false><<<nr_blocks, T, T * sizeof(float4)>>>(dt, N,
                (float4*) thrust::raw_pointer_cast(&ps_d[current_read]),
                (float4*) thrust::raw_pointer_cast(&ps_d[current_write]),
                NULL);
        cudaThreadSynchronize();

        std::swap(current_read, current_write);
    }

}

