#include "cuda-integrate.hpp"
#include <algorithm>

#define DEBUG 1
#if DEBUG
#    include "util/cuPrintf.cu"
#    include <iostream>
#    define MSG(s) (std::cout << s << std::endl)
#    define DBG(cmd) (cmd)
#else
#    define MSG(s) (void) 0
#    define DBG(cmd) (void) 0
#endif

#define PAD_TO 4

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

__constant__ float e2;
__constant__ float softening2;

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
    MSG("GPU arrays allocated (" << mem_size << " bytes)");

    cudaMemcpy(part_dev[0], &particles[0], mem_size, cudaMemcpyHostToDevice);
    MSG("Particle array copied to GPU");
    current_read = 0;
    current_write = 1;
}

void gpu_finalize() {
    cudaFree(vel_dev[0]);
    cudaFree(vel_dev[1]);
    cudaFree(part_dev[0]);
    cudaFree(part_dev[1]);
}

__device__ float gaussian_kernel_factor(float r2) {
    return (1 - expf(-r2/e2)) / (2 * M_PI * r2);
}


__device__ float2 vortex_interaction(float2 du, float4 self, float4 other) {

    float2 dif;
    dif.x = self.x - other.x;
    dif.y = self.y - other.y;

    float r2 = (dif.x * dif.x) + (dif.y * dif.y) + softening2;
    float kf = gaussian_kernel_factor(r2);
    du.x += -other.y * kf;
    du.y +=  other.x * kf;

    return du;
}


__device__ float2 tile_computation(float4 part, float2 du) {
    extern __shared__ float4 shared_part[];

    // Explanation for this on Nvidia's n-body program.
#ifdef _Win64
    unsigned long long i = 0;
#else
    unsigned long i = 0;
#endif

    for (unsigned int counter = 0; counter < blockDim.x; ) {
        du = vortex_interaction(du, part, shared_part[i++]);
        ++counter;
        // TODO: unroll loop
    }

    return du;
}


// divless mod
#define WRAP(x, m) ((x) < (m) ? (x) : (x) - (m))

__device__ float2 biot_savart_law(float4 part, float4* parts, unsigned nr_particles) {
    extern __shared__ float4 shared_part[];

    const unsigned p = blockDim.x;
    const unsigned b = blockIdx.x;
    const unsigned i = threadIdx.x;
    const unsigned nr_blocks = gridDim.x;

    int current_tile = b;
    float2 u = {0.0f, 0.0f};
#if 0
    for (int tile_count = 0; tile_count < nr_blocks; ++tile_count) {

        // fetch particle from global memory
        shared_part[i] = parts[p * b + i];
        __syncthreads();

        float2 du = {0.0f, 0.0f};
        du = tile_computation(part, du);
        u.x += du.x;
        u.y += du.y;
        __syncthreads();

        current_tile = WRAP(current_tile + 1, nr_blocks);
    }
#endif

    return u;
}


__global__ void
integrate(float4 *new_part, float4 *new_vel,
          float4 *old_part, float4 *old_vel, float dt, unsigned nr_particles) {
    unsigned index = blockIdx.x * blockDim.x + threadIdx.x;
    DBG(cuPrintf("Index: %u\n", index));

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
    size_t shared_mem_size = p * sizeof(float4);
    size_t grid_size = PAD(nr_particles) / p;

    MSG(nr_particles << " " << PAD(nr_particles) << " "  << p);
    DBG(cudaPrintfInit());
    for (int i = 0; i < nr_iterations; ++i) {
        MSG("integrate<<<" << grid_size << ", " << p << ", " << shared_mem_size << ">>>");
        integrate<<<grid_size, p, shared_mem_size>>>(
            (float4 *) part_dev[current_write], (float4 *) vel_dev[current_write],
            (float4 *) part_dev[current_read],  (float4 *) vel_dev[current_read],
            dt, PAD(nr_particles));

        std::swap(current_read, current_write);
    }
    DBG(cudaPrintfDisplay());
    DBG(cudaPrintfEnd());
}

