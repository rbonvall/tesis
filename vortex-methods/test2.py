#!/usr/bin/env python

from __future__ import division
from numpy import *

import init_position
import vm
import kernels
import problems.lamb_oseen

import functools
import sys

def p(s):
    print s
    sys.stdout.flush()

def main():
    import pylab

    x0, x1 = -0.5, 0.5
    y0, y1 = -0.5, 0.5
    h = .0625
    #h = .125
    nu = 5e-4
    dt = 0.01
    t0 = 0.01

    x, y = init_position.quasirandom(x0, x1, y0, y1, cell_size=h)

    # initial vorticity and circulation
    vort = problems.lamb_oseen.vorticity(x, y, t0, nu=nu)
    circ = h**2 * vort
    circ = vort/sum(vort)
    p("Total circulation: %f" % sum(circ))

    h_mesh = h
    M = int(ceil((x1 - x0)/h_mesh)) + 1
    N = int(ceil((y1 - y0)/h_mesh)) + 1
    blob_kernel = functools.partial(kernels.gaussian_cutoff, e2=h**2)
    i_mesh, j_mesh = mgrid[x0: x0 + M * h_mesh: h_mesh,
                           y0: y0 + N * h_mesh: h_mesh]
    w_mesh = vm.remesh_vorticity(x, y, circ, blob_kernel,
                                 x0, y0, h_mesh, M, N)

    plot_every = 100
    plot_rows, plot_cols = 2, 4

    t = t0
    iteration = 0
    plot_count = 0
    while True:
        plot_now = (iteration % plot_every == 0)

        eval_velocity = functools.partial(vm.eval_velocity, circ=circ, squared_blob_size=h**2)

        # begin Runge-Kutta integration
        u, v = zeros_like(x), zeros_like(y)

        kx, ky = eval_velocity(x, y) # k1
        u += kx/6
        v += ky/6

        dx, dy = kx * (dt/2), ky * (dt/2)
        kx, ky = eval_velocity(x + dx, y + dy) # k2
        u += kx/3
        v += ky/3

        dx, dy = kx * (dt/2), ky * (dt/2)
        kx, ky = eval_velocity(x + dx, y + dy) # k3
        u += kx/3
        v += ky/3

        dx, dy = kx * dt, ky * dt
        kx, ky = eval_velocity(x + dx, y + dy) # k4
        u += kx/6
        v += ky/6

        if plot_now:
            plot_count += 1
            pylab.subplot(plot_rows, plot_cols, plot_count)
            p("Starting plot %d" % plot_count)

            p("Remeshing vorticity at t = %f" % t)
            w_mesh = vm.remesh_vorticity(x, y, circ, blob_kernel,
                                         x0, y0, h_mesh, M, N)

            pylab.contour(i_mesh, j_mesh, w_mesh, 20)
            pylab.scatter(x, y, s=3)
            pylab.quiver(x, y, u, v)#, color='#444444', headwidth=2)

        x += u * dt
        y += v * dt

        iteration += 1
        t += dt

        #if plot_count == 2:
        if plot_count == plot_cols * plot_rows:
        #if plot_count == plot_cols:
            break

    pylab.show()


if __name__ == '__main__':
    main()

