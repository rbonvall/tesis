#!/usr/bin/env python
# vim: set fileencoding=utf-8:

from __future__ import division
from numpy import *

import init_position
import vm
import kernels
import vel_integration
import problems.lamb_oseen

import functools
import sys

def p(s):
    print s
    sys.stdout.flush()

def parse_args():
    import optparse

    parser = optparse.OptionParser()
    op = parser.add_option
    op('-x', type=float, nargs=2, default=(-0.5, 0.5), metavar='X₀ X₁')
    op('-y', type=float, nargs=2, default=(-0.5, 0.5), metavar='Y₀ Y₁')
    op('--cell-size', type=float, default=.0625, metavar='H')
    op('--viscosity', '-n', type=float, default=5e-4, metavar=u'ν')
    op('--t0', type=float, default=0.01, metavar=u'T₀')
    op('--dt', type=float, default=0.01, metavar=u'ΔT')
    op('--lattice',    action='store_const', dest='init', const=init_position.lattice,
                       help='Initialize particles in a lattice distribution',
                       default=init_position.lattice)
    op('--triangular', action='store_const', dest='init', const=init_position.triangular,
                       help='Initialize particles in a triangular distribution')
    op('--random',     action='store_const', dest='init', const=init_position.quasirandom,
                       help='Initialize particles in a quasirandom distribution')
    op('--circulation', type=float, default=1.0, metavar=u'Γ₀',
                       help='Lamb-Oseen total circulation')
    op('--plot-every', type=int, default=100, metavar='N',
                       help='Plot particles every N iteration')
    op('--euler', action='store_const', dest='timestepping', const=vel_integration.euler,
                        help='Use Euler time-stepping', default=vel_integration.euler)
    op('--runge-kutta', action='store_const', dest='timestepping', const=vel_integration.runge_kutta,
                        help='Use 4th order Runge-Kutta time-stepping')
    return parser.parse_args()


def main():
    import pylab

    options, args = parse_args()

    x0, x1 = options.x
    y0, y1 = options.y
    h = options.cell_size
    nu = options.viscosity
    t0 = options.t0
    dt = options.dt
    initialize = options.init
    plot_every = options.plot_every
    total_circulation = options.circulation
    velocity_integration = options.timestepping

    x, y = initialize(x0, x1, y0, y1, cell_size=h)

    # specialize vorticity functions
    vort_t = functools.partial(problems.lamb_oseen.vorticity, nu=nu)
    vort_0 = functools.partial(vort_t, t=t0)

    # initial vorticity and circulation
    vort = vort_0(x, y)
    circ = h**2 * vort
    circ = vort/sum(vort) * total_circulation

    h_mesh = h
    M = int(ceil((x1 - x0)/h_mesh)) + 1
    N = int(ceil((y1 - y0)/h_mesh)) + 1
    blob_kernel = functools.partial(kernels.gaussian_cutoff, e2=h**2)
    i_mesh, j_mesh = mgrid[x0: x0 + M * h_mesh: h_mesh,
                           y0: y0 + N * h_mesh: h_mesh]
    w_mesh = vm.remesh_vorticity(x, y, circ, blob_kernel,
                                 x0, y0, h_mesh, M, N)

    plot_rows, plot_cols = 2, 4

    t = t0
    iteration = 0
    plot_count = 0
    while True:
        plot_now = (iteration % plot_every == 0)

        u, v = velocity_integration(x, y, circ, dt, h**2)

        if plot_now:
            plot_count += 1
            pylab.subplot(plot_rows, plot_cols, plot_count)
            p("Starting plot %d" % plot_count)

            p("Remeshing vorticity at t = %f" % t)
            w_mesh = vm.remesh_vorticity(x, y, circ, blob_kernel,
                                         x0, y0, h_mesh, M, N)

            pylab.contour(i_mesh, j_mesh, w_mesh, 20)
            pylab.scatter(x, y, s=3)
            pylab.quiver(x, y, u, v)#, color='#444444', headwidth=2)

        x += u * dt
        y += v * dt

        iteration += 1
        t += dt

        #if plot_count == 2:
        if plot_count == plot_cols * plot_rows:
        #if plot_count == plot_cols:
            break

    pylab.show()


if __name__ == '__main__':
    main()

