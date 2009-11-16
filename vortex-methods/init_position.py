#!/usr/bin/env python

from __future__ import division
from numpy import *
import numpy.random

def lattice(x0, x1, y0, y1, cell_size):
    x, y = mgrid[x0 + cell_size/2 : x1 : cell_size,
                 y0 + cell_size/2 : y1 : cell_size]
    return x.flatten(), y.flatten()

def triangular(x0, x1, y0, y1, cell_size):
    x_even, y_even = mgrid[x0 + cell_size/2 : x1 : cell_size,
                           y0 + cell_size/2 : y1 : 2 * cell_size]
    x_odd, y_odd   = mgrid[x0 + cell_size     : x1 : cell_size,
                           y0 + 3*cell_size/2 : y1 : 2 * cell_size]
    x = hstack([x_even.flatten(), x_odd.flatten()])
    y = hstack([y_even.flatten(), y_odd.flatten()])
    return x, y

def quasirandom(x0, x1, y0, y1, cell_size):
    x, y = lattice(x0, x1, y0, y1, cell_size)
    x_noise = numpy.random.uniform(-cell_size/2, cell_size/2, size=x.shape)
    y_noise = numpy.random.uniform(-cell_size/2, cell_size/2, size=x.shape)
    return x + x_noise, y + y_noise


def main():
    import pylab

    x0, x1 = -1, 1
    y0, y1 = -1, 1
    cell_size = 0.25

    def new_subplot(n, x, y):
        s = pylab.subplot(130 + n, autoscale_on=False, xlim=(x0, x1), ylim=(y0, y1))
        #s.grid(True)
        pylab.scatter(x, y, s=1)

    x, y = lattice(x0, x1, y0, y1, cell_size)
    new_plot(1, x, y)

    x, y = triangular(x0, x1, y0, y1, cell_size)
    new_plot(2, x, y)

    x, y = quasirandom(x0, x1, y0, y1, cell_size)
    new_plot(3, x, y)

    pylab.show()

if __name__ == '__main__':
    main()
