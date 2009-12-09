#!/usr/bin/env python

from __future__ import division
import numpy
import numpy.random

def lattice(x0, x1, y0, y1, cell_size):
    x, y = numpy.mgrid[x0 + cell_size/2 : x1 : cell_size,
                       y0 + cell_size/2 : y1 : cell_size]
    return x.flatten(), y.flatten()

def triangular(x0, x1, y0, y1, cell_size):
    x_even, y_even = numpy.mgrid[x0 + cell_size/2 : x1 : cell_size,
                                 y0 + cell_size/2 : y1 : 2 * cell_size]
    x_odd, y_odd   = numpy.mgrid[x0 + cell_size     : x1 : cell_size,
                                 y0 + 3*cell_size/2 : y1 : 2 * cell_size]
    x = numpy.hstack([x_even.flatten(), x_odd.flatten()])
    y = numpy.hstack([y_even.flatten(), y_odd.flatten()])
    return x, y

def quasirandom(x0, x1, y0, y1, cell_size):
    x, y = lattice(x0, x1, y0, y1, cell_size)
    x_noise = numpy.random.uniform(-cell_size/2, cell_size/2, size=x.shape)
    y_noise = numpy.random.uniform(-cell_size/2, cell_size/2, size=x.shape)
    return x + x_noise, y + y_noise


def main():
    import pylab
    import optparse
    parser = optparse.OptionParser()
    parser.add_option('--x0', type=float, default=-1.0)
    parser.add_option('--x1', type=float, default=+1.0)
    parser.add_option('--y0', type=float, default=-1.0)
    parser.add_option('--y1', type=float, default=+1.0)
    parser.add_option('--cell-size', '-e', type=float, default=1e-1)
    (options, args) = parser.parse_args()

    x0, x1 = options.x0, options.x1
    y0, y1 = options.y0, options.y1
    cell_size = options.cell_size

    def new_subplot(n, x, y):
        s = pylab.subplot(130 + n, autoscale_on=False,
                          xlim=(x0, x1), ylim=(y0, y1))
        #s.grid(True)
        pylab.scatter(x, y, s=1)

    x, y = lattice(x0, x1, y0, y1, cell_size)
    new_subplot(1, x, y)

    x, y = triangular(x0, x1, y0, y1, cell_size)
    new_subplot(2, x, y)

    x, y = quasirandom(x0, x1, y0, y1, cell_size)
    new_subplot(3, x, y)

    pylab.show()

if __name__ == '__main__':
    main()
