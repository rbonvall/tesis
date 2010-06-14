#include "cuda-integrate.hpp"
#include <algorithm>

#define PAD_TO 256

#if defined(PAD_TO) && (PAD_TO > 0)
#    define PAD_B(x, b) ((x) % (b) == 0 ? (x) : (b) * (1 + (x) / (b)))
#    define PAD(x)      PAD_B((x), PAD_TO)
#else
#    define PAD(x)      (x)
#endif

// Two arrays for velocity and two for position
// will be allocated on the device:
// one of each for reading and one for writing.
// Variables current_{read,write} indicate which one is being used for what,
// and are swapped at the end of an iteration.
static float4 *vel_dev[2] = {NULL, NULL};
static float4 *part_dev[2] = {NULL, NULL};
static unsigned current_read;
static unsigned current_write;
static unsigned nr_particles;

void gpu_init(std::vector<particle>& particles) {
    nr_particles = particles.size();

    // convert particles to float4s for alignment purposes
    std::vector<p4> particles4;
    particles4.reserve(PAD(nr_particles));

    for (int i = 0; i < nr_particles; ++i) {
        particles4.push_back(p4(particles[i]));
    }
    for (int i = nr_particles; i < PAD(nr_particles); ++i) {
        particles4.push_back(p4(0.0, 0.0, 0.0));
    }

    unsigned mem_size = sizeof(float) * 4 * PAD(nr_particles);
    cudaMalloc((void **) &part_dev[0], mem_size);
    cudaMalloc((void **) &part_dev[1], mem_size);
    cudaMalloc((void **) &vel_dev[0], mem_size);
    cudaMalloc((void **) &vel_dev[1], mem_size);

    cudaMemcpy(part_dev[0], &particles[0], mem_size, cudaMemcpyHostToDevice);
    current_read = 0;
    current_write = 1;
}

void gpu_finalize() {
    cudaFree(vel_dev[0]);
    cudaFree(vel_dev[1]);
    cudaFree(part_dev[0]);
    cudaFree(part_dev[1]);
}


__device__ float2 biot_savart_law(float4 p, float4* pos, int nr_particles) {
    extern __shared__ float4 shared_pos[];

    // TODO: implement B-S law
    float2 u;
    u.x = 1; u.y = 1;
    return u;
}


__global__ void
integrate(float4 *new_part, float4 *new_vel,
          float4 *old_part, float4 *old_vel, float dt, unsigned nr_particles) {
    unsigned index = blockIdx.x * blockDim.x + threadIdx.x;

    // copy particle from global memory
    float4 p = old_part[index];

    // compute velocity by applying B-S law along the tile
    float2 vel = biot_savart_law(p, old_part, nr_particles);

    // convect particle and copy it to global memory
    p.x += vel.x * dt;
    p.y += vel.y * dt;
    new_part[index] = p;

    // copy computed velocity to global memory
    new_vel[index].x = vel.x;
    new_vel[index].y = vel.y;
}

void vm_integrate(float dt, unsigned nr_iterations, int p) {
    int shared_mem_size = p * sizeof(float4);
    for (int i = 0; i < nr_iterations; ++i) {
        integrate<<<nr_particles / p, p, shared_mem_size>>>(
            (float4 *) part_dev[current_write], (float4 *) vel_dev[current_write],
            (float4 *) part_dev[current_read],  (float4 *) vel_dev[current_read],
            dt, PAD(nr_particles));

        std::swap(current_read, current_write);
    }
}

