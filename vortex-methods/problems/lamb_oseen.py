#!/usr/bin/env python

from __future__ import division
from numpy import *

def vorticity(x, y, t, gamma0=1e0, nu=5e-4):
    r_sq = x**2 + y**2
    den = (4 * nu * t)
    return (gamma0/(den * pi)) * exp(-r_sq/den)

def velocity(x, y, t, gamma0=1e0, nu=5e-4):
    r_sq = x**2 + y**2
    r = hypot(x, y)
    c = gamma0/(2 * pi * r)
    kf = c * (1 - exp(-r_sq/(4 * nu * t)))
    return -y * kf, x * kf


def main():
    import pylab

    N = 2**8
    h = 1/N
    dom = linspace(-0.1, 0.1, N + 1)
    x, y = meshgrid(dom, dom)
    x_vel, y_vel = x[::16, ::16], y[::16, ::16]

    plot_rows, plot_cols = 2, 5

    #w0 = h**2 * lamb_oseen(x, y, 0.01)  # VM:T&P, page 27

    t0 = 0.1
    dt = 0.1

    t = t0
    iteration = 0
    for iteration in range(plot_cols - 1):

        pylab.subplot(plot_rows, plot_cols, iteration + 1)
        w = vorticity(x, y, t)
        if iteration == 0:
            cs = pylab.contour(x, y, w)
            l = cs.levels
        else:
            cs = pylab.contour(x, y, w, l)
        pylab.colorbar()

        pylab.subplot(plot_rows, plot_cols, iteration + 1 + plot_cols)
        u, v = velocity(x_vel, y_vel, t)
        vel = hypot(u, v)
        pylab.quiver(x_vel, y_vel, u, v)
        pylab.contour(x_vel, y_vel, vel)


        iteration += 1
        t += dt

        pylab.subplot(plot_rows, plot_cols, plot_rows * plot_cols)
        pylab.plot(dom, vorticity(dom, 0, t))
        pylab.title('Vorticity decay')

    pylab.show()


if __name__ == '__main__':
    main()
