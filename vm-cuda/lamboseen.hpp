#pragma once
#include <cmath>

struct lamb_oseen_vortex {
    float total_circulation;  // \gamma_0
    float viscosity;          // \nu

    lamb_oseen_vortex(float circ, float visc) :
            total_circulation(circ), viscosity(visc) {}

    float operator() (float x, float y, float t) {
        float one_over_four_nu_t = 1/(4 * t * viscosity);
        float r_squared = x * x + y * y;
        return total_circulation * one_over_four_nu_t *
               exp(-r_squared * one_over_four_nu_t) / M_PI;
    }

    float velocity(float x, float y, float t) {
        float r_squared = x * x + y * y;
        return total_circulation *
               (1 - exp(-r_squared / sqrt(4 * t * viscosity))) /
               (2 * M_PI * sqrt(r_squared));
    }
};

