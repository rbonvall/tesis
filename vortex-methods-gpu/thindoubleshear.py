#!/usr/bin/env python

from __future__ import division
from vm import *
from pylab import *

def u0(x, y, rho=80, delta=0.05):
    return tanh(rho * minimum(y - 0.25, 0.75 - y))

def v0(x, y, rho=80, delta=0.05):
    return delta * sin(2 * pi * (x + 0.25))

dt = 2e-2
t_end = 10e-2
M = 64
h = 1/M
viscosity = 1e-4

x_dom = linspace(0.0, 1.0, M)
y_dom = linspace(0.0, 1.0, M)
x_m, y_m = meshgrid(x_dom, y_dom)

u_m, v_m = u0(x_m, y_m), v0(x_m, y_m)
w_m = vorticity(u_m, v_m)

figure(1)
s = (slice(None, None, M//32),) * 2
for t in arange(0.0, t_end, dt):

    subplot(341)
    contourf(x_m, y_m, w_m)
    title('Initial vorticity')

    # integrate vorticity
    print "DIFFUSION"
    w_m += dt * viscosity * diffusion(w_m, h)
    print "RESHAPE"
    w_p = w_m.reshape(w_m.size)

    subplot(345)
    contourf(x_m, y_m, w_m)
    title('Vorticity after diffusion')

    subplot(342)
    contourf(x_m, y_m, u_m)
    title('Initial u')

    subplot(346)
    contourf(x_m, y_m, v_m)
    title('Initial v')

    subplot(3,4,10)
    quiver(x_m[s], y_m[s], u_m[s], v_m[s])
    title('Initial (u, v)')

    # integrate velocity
    print "POISSON SOURCE"
    f, g = poisson_source(w_m)

    subplot(343)
    contourf(x_m, y_m, f)
    title('Poisson source f')

    subplot(347)
    contourf(x_m, y_m, g)
    title('Poisson source g')

    print "POISSON SOLVE"
    u_m, v_m = poisson_solve(f, g, h)

    subplot(344)
    contourf(x_m, y_m, u_m)
    title('Solved u')

    subplot(348)
    contourf(x_m, y_m, v_m)
    title('Solved v')

    subplot(3,4,12)
    quiver(x_m[s], y_m[s], u_m[s], v_m[s])
    title('Solved (u, v)')

    # final particles' positions
    print "INTEGRATE"
    x_p = (x_m + u_m * dt).reshape(x_m.size)
    y_p = (y_m + v_m * dt).reshape(y_m.size)

    # remesh vorticity
    #print "REMESH"
    #w_m = remesh(w_p, x_p, y_p, x_m, y_m, h)

    subplot(349)
    contourf(x_m, y_m, w_m)
    title('Remeshed vorticity')

    show()
    break

