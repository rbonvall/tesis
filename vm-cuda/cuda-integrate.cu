#include "cuda-integrate.hpp"

// Two arrays for velocity and two for position
// will be allocated on the device:
// one of each for reading and one for writing.
// Variables current_{read,write} indicate which one is being used for what,
// and are swapped at the end of an iteration.
static float *vel_dev[2] = {NULL, NULL};
static float *pos_dev[2] = {NULL, NULL};
static unsigned current_read = 0;
static unsigned current_write = 1;

void gpu_init(unsigned nr_particles) {
    unsigned mem_size = sizeof(float) * 4 * nr_particles;
    cudaMalloc((void **) vel_dev[0], mem_size);
    cudaMalloc((void **) vel_dev[1], mem_size);
    cudaMalloc((void **) pos_dev[0], mem_size);
    cudaMalloc((void **) pos_dev[1], mem_size);

}

void gpu_finalize() {
    cudaFree(vel_dev[0]);
    cudaFree(vel_dev[1]);
    cudaFree(pos_dev[0]);
    cudaFree(pos_dev[1]);
}

/* CUDA kernels */

__global__ void
integrate() {

}

void vm_integrate(std::vector<float>& x, std::vector<float>& y,
                  std::vector<float>& circ,
                  std::vector<float>& u, std::vector<float>& v) {

    unsigned nr_particles = x.size();

    /* copy from host to device */

//    integrate<<< >>>(...);


    /* copy from device to host */
}

