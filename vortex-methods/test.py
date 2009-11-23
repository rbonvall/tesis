#!/usr/bin/env python

from __future__ import division
from numpy import *
import pylab

import init_position
import vm
import problems.lamb_oseen


def main():
    x0, x1 = -0.5, 0.5
    y0, y1 = -0.5, 0.5
    h = .125
    nu = 1.0
    dt = 0.1
    t0 = 0.01

    plot_rows, plot_cols = 3, 3
    plot_every = 10

    x, y = init_position.quasirandom(x0, x1, y0, y1, cell_size=h)

    # initial vorticity and circulation
    vort = problems.lamb_oseen.vorticity(x, y, t0, nu=nu)
    circ = h**2 * vort

    # particle colors
    from itertools import cycle, izip as zip
    colors = 'red green blue magenta cyan yellow'.split()
    colors = [c for (c, _) in zip(cycle(colors), x)]

    t = t0
    iteration = 0
    plot_count = 0
    while True:
        plot_now = (iteration % plot_every == 0)

        if plot_now:
            plot_count += 1
            pylab.subplot(plot_rows, plot_cols, plot_count,
                          autoscale_on=False,
                          xlim=(1.2 * x0, 1.2 * x1),
                          ylim=(1.2 * y0, 1.2 * y1))
            pylab.scatter(x, y, s=2, c=colors, edgecolors='none')
            pylab.grid(True)

        u, v = vm.eval_velocity(x, y, circ)

        #pylab.quiver(x[::10], y[::10], u[::10], v[::10])

        if plot_now:
            pylab.quiver(x, y, u, v, color=colors)
            pylab.title('$t = %.2f$' % t)

        # convect
        x += u * dt
        y += v * dt

        iteration += 1
        t += dt

        if plot_count == plot_rows * plot_cols:
            break

    pylab.suptitle(ur'Lamb Oseen vortex, $\nu = %.3f$' % nu)
    pylab.show()


if __name__ == '__main__':
    main()
