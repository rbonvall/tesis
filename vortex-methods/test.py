#!/usr/bin/env python

from __future__ import division
from numpy import *
import pylab

import init_position
import vm
import lamb_oseen


def main():
    x0, x1 = -1, 1
    y0, y1 = -1, 1
    h = .25
    nu = 1e3
    t_f = 1.0
    dt = 1e-2

    plot_grid_size = 3, 3
    plot_every = 10

    x, y = init_position.lattice(x0, x1, y0, y1, cell_size=h)

    # initial vorticity and circulation
    vort = lamb_oseen.lamb_oseen(x, y, 0, nu=nu)
    circ = h**2 * w

    for it, t in enumerate(arange(t_f + dt, step=dt)):
        u, v = vm.eval_velocity(x, y, circ)
        x += u * dt
        y += v * dt

if __name__ == '__main__':
    main()
