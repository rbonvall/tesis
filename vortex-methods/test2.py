#!/usr/bin/env python

from __future__ import division
from numpy import *

import init_position
import vm
import kernels
import problems.lamb_oseen

import functools


def main():
    import pylab

    x0, x1 = -0.5, 0.5
    y0, y1 = -0.5, 0.5
    h = .0625
    nu = 1e3
    dt = 0.01
    t0 = 0.01

    plot_rows, plot_cols = 3, 3
    plot_every = 10

    x, y = init_position.triangular(x0, x1, y0, y1, cell_size=h)

    # initial vorticity and circulation
    vort = problems.lamb_oseen.vorticity(x, y, t0, nu=nu)
    circ = h**2 * vort

    blob_kernel = functools.partial(kernels.gaussian_cutoff, e2=0.1)
    w_mesh = vm.remesh_vorticity(x, y, circ, blob_kernel,
                                 x0=-1.0, y0=-1.0, h=0.1, M=20, N=20)

    pylab.contour(w_mesh)




if __name__ == '__main__':
    main()

