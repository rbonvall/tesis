#pragma once

#include <cuda_runtime_api.h>
#include <vector>

void vm_integrate(std::vector<float>& x, std::vector<float>& y,
                  std::vector<float>& circ,
                  std::vector<float>& u, std::vector<float>& v);

void gpu_init(unsigned nr_particles);
void gpu_finalize();
