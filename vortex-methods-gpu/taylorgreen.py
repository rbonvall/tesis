#!/usr/bin/env python

from __future__ import division
from numpy import sin, cos, pi, exp, arange, meshgrid, gradient
import matplotlib
matplotlib.use('Agg')
from pylab import *
from functools import partial
from optparse import OptionParser
from os import system


def u(x, y, t, U=0.5, Re=1e6):
    b = -8 * (pi**2) / Re
    return -U * exp(b * t) * cos(2 * pi * x) * sin(2 * pi * y)

def v(x, y, t, U=0.5, Re=1e6):
    b = -8 * (pi**2) / Re
    return  U * exp(b * t) * sin(2 * pi * x) * cos(2 * pi * y)

def main(t0, dt, t_end, dx, output_filename):
    dom = arange(0, 1, 0.1)
    dom = arange(0, 1, 0.05)
    x, y = meshgrid(dom, dom)

    #t0, dt = 0.0, 1e6
    #t0, dt, t_end = 0.0, 2e-4, 0.1
    t0, dt, t_end = 0.0, 10.0, 2000
    i, t = 0, t0
    while t < t_end:
        figure(i)
        t = t0 + i * dt
        ut, vt = u(x, y, t), v(x, y, t)
        ux, uy = gradient(ut)
        vx, vy = gradient(vt)
        rot = uy - vx
        contour(x, y, rot)
        quiver(x, y, ut, vt)
        title(r't = %g' % t)
        savefig('tg-%05d.png' % i)
        if i % 10 == 0:
            print i
        i += 1
    system("mencoder 'mf://tg-?????.png' -mf type=png:fps=10 -ovc lavc "
           "-lavcopts vcodec=wmv2 -oac copy -o %s" % output_filename)
    system("rm tg-?????.png")

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('--t0', type="float", default=0.0)
    parser.add_option('--dt', type="float", default=0.1)
    parser.add_option('--t-end', type="float", default=1.0)
    parser.add_option('--dx', type="float", default=1.0)
    parser.add_option('-o', type="string", default='tg.avi')
    (options, args) = parser.parse_args()

    main(t0=options.t0,
         dt=options.dt,
         t_end=options.t_end,
         dx=options.dx,
         output_filename=options.o)
