#!/usr/bin/env python
# encoding: utf-8

from __future__ import division
from numpy import *
from scipy.linalg import norm
import init_position
import pylab

class LambOseenVortex:
    '''
    Analytical description of the Lamb-Oseen vortex.

    >>> vortex = LambOseenVortex()
    >>> w = vortex(2, 3, 1)             # vorticity in x = (2, 3) at t = 1
    >>> u, v = vortex.velocity(2, 3, 1) # vorticity in x = (2, 3) at t = 1
    '''

    def __init__(self, total_circulation=1.0, viscosity=5e-4):
        self.total_circulation = total_circulation
        self.viscosity = viscosity

    def __call__(self, x, y, t):
        one_over_four_nu_t = 1 / (4 * t * self.viscosity)
        r_squared = x * x + y * y
        return (self.total_circulation * one_over_four_nu_t *
                exp(-r_squared * one_over_four_nu_t) / pi)

    def velocity(self, x, y, t):
        r = hypot(x, y)
        r_squared = r ** 2
        factor = (self.total_circulation *
                  (1 - exp(-r_squared / (4 * t * self.viscosity))) /
                  (2 * pi * r))
        return (-y * factor, x * factor)

# ζ_ε(r²)
def gaussian_cutoff(r2, e2):
    c = 1/(pi * e2)
    return c * exp(-r2/e2)

# K_ε(r²)
def gaussian_bs_kernel_factor(r2, e2):
    return 1 / (2 * pi * r2) * (1 - exp(-r2/e2))

# K(r²)
def bs_kernel_factor(r2):
    return 1 / (2 * pi * r2)

# η_ε(r²)
diffusion_kernel = gaussian_cutoff


# Biot-Savart law
def eval_velocity(x, y, circ, e2):
    u, v = zeros_like(x), zeros_like(y)
    for p, (x_p, y_p) in enumerate(zip(x, y)):
        dx, dy = x_p - x, y_p - y
        r2 = dx**2 + dy**2
        kf = gaussian_bs_kernel_factor(r2, e2)
        K1, K2 = -kf * dy, kf * dx
        K1[p], K2[p] = 0.0, 0.0
        u[p], v[p] = dot(circ, K1), dot(circ, K2)
    return u, v

# Particle strength exchange
def eval_circulation_change(x, y, circ, nu, e2):
    dcirc = zeros_like(circ)
    for p, (x_p, y_p, circ_p) in enumerate(zip(x, y, circ)):
        dx, dy = x - x_p, y - y_p
        r2 = dx**2 + dy**2
        eta = diffusion_kernel(r2, e2)
        circ_diff = circ - circ_p
        circ[p] = (nu * (e2/4) / e2) * dot(circ_diff, eta)

        # particle volume has been approximated
        # as v_p = h² ≈ ε²/4
    return dcirc


# compares analytical and computed velocity at particle locations
def squared_velocity_errors(x, y, u, v, t, vortex):
    # analytical velocity
    u_a, v_a = vortex.velocity(x, y, t)
    return ((u - u_a) ** 2 + (v - v_a) ** 2)

def main():
    nu = 5e-4     # ν
    gamma0 = 1.0  # Γ₀
    vortex = LambOseenVortex(total_circulation=gamma0, viscosity=nu)

    # particle initialization
    x0, x1 = -0.5, 0.5
    y0, y1 = -0.5, 0.5
    h = .125 / 4
    e = 2 * h     # ε
    x, y = init_position.triangular(x0, x1, y0, y1, cell_size=h)

    # integration parameters
    t0 = 4
    dt = 0.001      # δt
    nr_iterations = 100

    # initial circulation evaluation
    initial_vorticity = vortex(x, y, t0)
    circ = h**2 * initial_vorticity

    iteration = 0
    t = t0
    while True:

        # convection
        u, v = eval_velocity(x, y, circ, e ** 2)
        x += u * dt
        y += v * dt

        # diffusion
        #dcirc = eval_circulation_change(x, y, circ, nu, e ** 2)
        #circ += dcirc * dt

        if iteration % 10 == 0:
            err2 = squared_velocity_errors(x, y, u, v, t, vortex)
            print "Iteration {0}, t = {1}".format(iteration, t)
            print "Squared error: {0} ± {1}".format(err2.mean(), err2.std())
            print

            #pylab.subplot(2, 5, (iteration + 1) // 10)
            #u_a, v_a = vortex.velocity(x, y, t)
            #pylab.scatter(x, y, s=1, c='red', edgecolor='none')
            #pylab.quiver(x, y, u, v, color='red')
            #pylab.quiver(x, y, u_a, v_a, color='blue')

        if not (iteration < nr_iterations):
            break

        iteration += 1
        t += dt

    u_a, v_a = vortex.velocity(x, y, t)
    pylab.scatter(x, y, s=1, c='red', edgecolor='none')
    pylab.quiver(x, y, u, v, color='red')
    pylab.quiver(x, y, u_a, v_a, color='blue')

    pylab.show()


def test_velocity_evaluation():
    x    = array([-2e-2,  1e-2, 2e-2])
    y    = array([ 1e-2, -1e-2, 2e-2])
    circ = array([ 1e-3,  1e-3, 2e-3])
    e = 3e-2

    # expected results with precision 1e-5
    u   = array([-0.28e-3, 8.28e-3, -4.00e-3])
    v   = array([-9.16e-3, 0.66e-3,  4.25e-3])

    # obtained results
    u_o, v_o = eval_velocity(x, y, circ, e**2)

    if norm(hypot(u_o - u, v_o - v)) < 5e-5:
        print 'Test of velocity evaluation failed:'
        print 'u expected: ', u
        print 'v expected: ', v
        print 'u obtained: ', u_o
        print 'v obtained: ', v_o

if __name__ == '__main__':
    test_velocity_evaluation()
    main()

