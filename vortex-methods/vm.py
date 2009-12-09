#!/usr/bin/env python

from __future__ import division
from numpy import *
from numpy.linalg import norm
from kernels import gaussian_bs_kernel_factor as gaussian

def eval_velocity(x, y, circ, squared_blob_size=1.0, bs_kernel_factor=gaussian):
    u, v = zeros_like(x), zeros_like(y)
    for p, (x_p, y_p) in enumerate(zip(x, y)):
        r2 = (x - x_p)**2 + (y - y_p)**2
        kf = bs_kernel_factor(r2, squared_blob_size)
        K1, K2 = -kf * y, kf * x
        K1[p], K2[p] = 0.0, 0.0
        u[p], v[p] = dot(circ, K1), dot(circ, K2)
    return u, v


def remesh_vorticity(x, y, circ, blob_kernel, x0, y0, h, M, N):

    s = 16 # TODO: compute in terms of the blob kernel support
    w_mesh = zeros((M, N))

    for p, (x_p, y_p, circ_p) in enumerate(zip(x, y, circ)):
        i_p, j_p = int(floor((x_p - x0)/h)), int(floor((y_p - y0)/h))

        # create computational support
        i_supp_start, i_supp_end = max(0, i_p - s//2 + 1), min(i_p + s//2, M - 1) + 1
        j_supp_start, j_supp_end = max(0, j_p - s//2 + 1), min(j_p + s//2, N - 1) + 1
        i_supp, j_supp = mgrid[i_supp_start:i_supp_end, j_supp_start:j_supp_end]
        x_supp, y_supp = x0 + i_supp * h, y0 + j_supp * h

        # evaluate vorticity
        dx, dy = x_supp - x_p, y_supp - y_p
        weights = blob_kernel(dx**2 + dy**2)
        w_mesh[i_supp, j_supp] += circ_p * weights

    return w_mesh


def main():
    pass

if __name__ == '__main__':
    main()
