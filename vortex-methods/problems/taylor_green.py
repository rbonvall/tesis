#!/usr/bin/env python

from __future__ import division
from numpy import *

def velocity(x, y, t, U=0.5, Re=1e6):
    nu = 1 / Re
    b = -8 * nu * pi**2

    f = U * exp(b * t)
    tpx, tpy = 2 * pi * x, 2 * pi * y
    return (-f * cos(tpx) * sin(tpy),
             f * sin(tpx) * cos(tpy))

def vorticity(x, y, t, U=0.5, Re=1e6):
    nu = 1 / Re
    b = -8 * nu * pi**2
    
    f = U * exp(b * t) * 4 * pi
    return f * cos(2 * pi * x, 2 * pi * y)


def main():
    import pylab

    N = 2**10
    h = 1/N
    dom = linspace(0.0, 1.0, N + 1)
    x, y = meshgrid(dom, dom)
    x_vel, y_vel = x[::16, ::16], y[::16, ::16]

    plot_rows, plot_cols = 2, 5

    t0 = 0.1
    dt = 0.5

    t = t0
    iteration = 0
    for iteration in range(plot_cols - 1):

        pylab.subplot(plot_rows, plot_cols, iteration + 1)
        values = vorticity(x, y, t)
        #if iteration == 0:
        #    cs = pylab.contour(x, y, values)
        #    l = cs.levels
        #else:
        #    cs = pylab.contour(x, y, values, l)
        pylab.contour(x, y, values)
        pylab.colorbar()

        pylab.subplot(plot_rows, plot_cols, iteration + 1 + plot_cols)
        u, v = velocity(x_vel, y_vel, t)
        vel = hypot(u, v)
        pylab.quiver(x_vel, y_vel, u, v)
        pylab.contour(x_vel, y_vel, vel)


        iteration += 1
        t += dt

        pylab.subplot(plot_rows, plot_cols, plot_rows * plot_cols)
        pylab.plot(dom, vorticity(dom, zeros_like(dom), t))
        pylab.title('Vorticity decay')

    pylab.show()

if __name__ == '__main__':
    main()
