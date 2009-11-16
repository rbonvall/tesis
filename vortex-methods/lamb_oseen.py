#!/usr/bin/env python

from __future__ import division
from numpy import *
import pylab

def lamb_oseen(x, y, t, gamma0=1e0, nu=1e3):
    r_sq = x**2 + y**2
    den = (4 * nu * t)
    return (gamma0/(den * pi)) * exp(-r_sq/den)


def plot_lamb_oseen():
    N = 2**8
    h = 1/N
    dom = linspace(-1, 1, N + 1)
    x, y = meshgrid(dom, dom)
    #
    w0 = h**2 * lamb_oseen(x, y, 0.01)  # VM:T&P, page 27
    #
    for i in range(6):
        pylab.subplot(2, 6, i + 1)
        t = .01 + (i + 1) * .0001
        if i == 0:
            cs = pylab.contour(x, y, lamb_oseen(x, y, t))
            l = cs.levels
        else:
            cs = pylab.contour(x, y, lamb_oseen(x, y, t), l)
        pylab.colorbar()
        pylab.subplot(2, 6, i + 7)
        pylab.plot(dom, lamb_oseen(dom, 0, t))
    pylab.show()

def main():
    pass


if __name__ == '__main__':
    plot_lamb_oseen()
    #main()
