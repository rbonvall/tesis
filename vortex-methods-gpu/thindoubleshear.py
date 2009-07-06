#!/usr/bin/env python

from __future__ import division
from vm import *

def u0(x, y, rho=80, delta=0.05):
    return tanh(rho * minimum(y - 0.25, 0.75 - y))

def v0(x, y, rho=80, delta=0.05):
    return delta * sin(2 * pi * (x + 0.25))

dt = 2e-2
t_end = 10e-2
M = 64
h = 1/M

x_dom = linspace(0.0, 1.0, M + 1)
y_dom = linspace(0.0, 1.0, M + 1)
x_m, y_m = meshgrid(x_dom, y_dom)

u_m, v_m = u0(x_m, y_m), v0(x_m, y_m)
w_m = vorticity(u_m, v_m)

for t in arange(0.0, t_end, dt):
    # integrate vorticity
    print "DIFFUSION"
    w_m += diffusion(w_m, h)
    print "RESHAPE"
    w_p = w_m.reshape(w_m.size)

    # integrate velocity
    print "POISSON SOURCE"
    f, g = poisson_source(w_m)
    print "POISSON SOLVE"
    u, v = poisson_solve(f, g, h)

    # final particles' positions
    print "INTEGRATE"
    x_p = (x_m + u_m * dt).reshape(x_m.size)
    y_p = (y_m + v_m * dt).reshape(y_m.size)

    # remesh vorticity
    print "REMESH"
    w_mesh = remesh(w_p, x_p, y_p, x_m, y_m, h)

    print t,

    


contourf(x_m, y_m, rot)
quiver(x_m, y_m, u0, v0)

show()


        


