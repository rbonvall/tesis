#pragma once

#include <cuda_runtime_api.h>
#include <vector>
#include "vm.hpp"

// same structure of a GPU float4
struct p4 {
    float x, y, circ, padding;

    p4(float x, float y, float circ) :
        x(x), y(y), circ(circ), padding(0.0) {}

    p4(particle& p) :
        x(p.x), y(p.y), circ(p.circ), padding(0.0) {}
};

void vm_integrate(std::vector<float>& x, std::vector<float>& y,
                  std::vector<float>& circ,
                  std::vector<float>& u, std::vector<float>& v);

void vm_integrate(float dt, unsigned nr_iterations = 1, int p = 256);
void gpu_init(std::vector<particle>& particles);

void gpu_init(unsigned nr_particles);
void gpu_finalize();
