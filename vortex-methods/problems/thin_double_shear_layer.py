#!/usr/bin/env python

from __future__ import division
from numpy import *

def velocity(x, y, delta=0.05, rho=80, nu=1e-4):
    u = tanh(rho * minimum(y - 0.25, 0.75 - y))
    v = delta * sin(2 * pi * (x + 0.25))
    return u, v


def main():
    import pylab

    N = 2**4
    h = 1/N
    x, y = mgrid[0.0:1.0:h, 0.0:1.0:h]

    pylab.quiver(x, y, *velocity(x, y))
    pylab.show()


if __name__ == '__main__':
    main()
