#!/usr/bin/env python

from __future__ import division
from numpy import *
import pylab

def lamb_oseen(x, y, t, gamma0=1e0, nu=5e-4):
    r_sq = x**2 + y**2
    den = (4 * nu * t)
    return (gamma0/(den * pi)) * exp(-r_sq/den)


def plot_lamb_oseen():
    N = 2**8
    h = 1/N
    dom = linspace(-0.1, 0.1, N + 1)
    x, y = meshgrid(dom, dom)

    plot_rows, plot_cols = 2, 4

    #w0 = h**2 * lamb_oseen(x, y, 0.01)  # VM:T&P, page 27

    t0 = 0.1
    dt = 0.1

    t = t0
    iteration = 0
    for iteration in range(4):

        pylab.subplot(plot_rows, plot_cols, iteration + 1)
        lamb_oseen_values = lamb_oseen(x, y, t)
        if iteration == 0:
            cs = pylab.contour(x, y, lamb_oseen_values)
            l = cs.levels
        else:
            cs = pylab.contour(x, y, lamb_oseen_values, l)
        pylab.colorbar()

        pylab.subplot(plot_rows, plot_cols, iteration + 1 + plot_cols)
        pylab.plot(dom, lamb_oseen(dom, 0, t))

        iteration += 1
        t += dt

    pylab.show()

def main():
    pass


if __name__ == '__main__':
    plot_lamb_oseen()
    #main()
