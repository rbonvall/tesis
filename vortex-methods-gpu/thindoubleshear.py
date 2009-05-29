#!/usr/bin/env python

from __future__ import division
from vm import *

def u0(x, y, rho=80, delta=0.05):
    return tanh(rho * minimum(y - 0.25, 0.75 - y))

def v0(x, y, rho=80, delta=0.05):
    return delta * sin(2 * pi * (x + 0.25))

dt = 2e-2
t_end = 1.0
M = 20
h = 1/M

# initial conditions
dom = linspace(0, 1, M + 1)
x_mesh, y_mesh = meshgrid(dom, dom)

# initial particle positions
x, y = x_mesh.copy().flatten(), y_mesh.copy().flatten()

for t in arange(0.0, t_end, dt):
    u, v = u0(x, y), v0(x, y)
    w = vorticity(u, v)

    # diffusion?

    # final particles' positions
    x = x + u * dt
    y = y + v * dt

    # remesh vorticiy
    w_mesh = remesh(w, x, y, x_mesh, y_mesh, h)

    


    




    





contourf(x, y, rot)
quiver(x, y, u0, v0)

show()


        


