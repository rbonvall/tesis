#!/usr/bin/env python

from __future__ import division
from numpy import pi, piecewise, exp

def exact_bs_kernel_factor(r2):
    return 1/(2 * pi * r2)


# Radial-symmetric cutoffs as a function of
# distance to origin and cutoff radius.

# Kernel factors are the factor of (-y, x) in
# the curl of the Green's function of the Laplacian.

def tophat_cutoff(r2, e2):
    return piecewise(r2, [r2 <= e2], [1/(pi*e2)])

def tophat_bs_kernel_factor(r2, e2):
    #e2 = ones_like(r2) * e2
    return 1/(2 * pi * piecewise(r2, [r2 <= e2, r2 > e2],
                                     [lambda x: e2, lambda x: x]))


def p1_cutoff(r2, e2):
    c = 2/pi
    return c * (2 - r2) / (1 - r2)**4

def p1_bs_kernel_factor(r2, e2):
    c = 1 / (2*pi)
    num = r2**2 + 3 * e2 * r2 + 4 * e2**2
    den = (e2 + r2)**3
    return c * num/den


def gaussian_cutoff(r2, e2):
    c = 1/(2 * pi * e2)
    return c * exp(-r2/(2 * e2))

def gaussian_bs_kernel_factor(r2, e2):
    c = 1/(2 * pi * r2)
    return c * (1 - exp(-r2/(2 * e2)))


def p2_e_cutoff(r2, e2):
    c = 1/(2 * pi)
    return c * (4 - r2) * exp(-r2)

def p2_e_bs_kernel_factor(r2, e2):
    c = 1/(2 * pi * r2)
    f = r2/e2
    return c * (1 + (f - 1) * exp(-f))


def p4_e_cutoff(r2, e2):
    c = 1/pi
    return (6 - 6 * r2 + r2**2) * exp(-r2)

def p4_e_bs_kernel_factor(r2, e2):
    c = 1/(2 * pi * r2)
    f = r2/e2
    return c * (1 - (1 - 2 * f + f**2/2) * exp(-f))



def main():
    import numpy
    import pylab

    e2 = 1.0

    # plot layout
    M, N = 2, 4

    # parameters for cutoff plotting
    r = numpy.linspace(0.0, 3.0, 100)
    r2 = r**2

    # parameters for Biot-Savart kernel plotting
    X = 5.0
    x_grid, y_grid = numpy.mgrid[-X:X:0.5,
                                 -X:X:0.5]
    r2_grid = x_grid**2 + y_grid**2


    def sp(n):
        pylab.subplot(M, N, n)
                      #autoscale_on=False,
                      #xlim=(r[0], r[-1]),
                      #ylim=(-0.05, .05+1/(pi*e2)))

    def cutoff_plot(f):
        pylab.plot(r, f(r2, e2))
        pylab.title(f.func_name)

    def bs_kernel_plot(f):
        kf = f(r2_grid, e2)
        pylab.quiver(x_grid, y_grid, -kf * y_grid, kf * x_grid)
        pylab.title(f.func_name.replace('_factor', ''))

    fns = [
        (tophat_cutoff,   tophat_bs_kernel_factor),
        (gaussian_cutoff, gaussian_bs_kernel_factor),
        (p2_e_cutoff,     p2_e_bs_kernel_factor),
        (p4_e_cutoff,     p4_e_bs_kernel_factor),
    ]

    for n, (cutoff, bs_kernel) in enumerate(fns):
        sp(n + 1);
        cutoff_plot(cutoff)
        sp(N + n + 1);
        bs_kernel_plot(bs_kernel)

    pylab.show()

if __name__ == '__main__':
    main()
