#!/usr/bin/env python

from __future__ import division
from numpy import *
from numpy.linalg import norm
from kernels import gaussian_bs_kernel_factor as gaussian

def biot_savart_kernel(x, y):
    denominator = 2 * pi * (x**2 + y**2)
    return (-y/denominator, x/denominator)

def eval_velocity(x, y, circ, squared_blob_size=1.0, bs_kernel_factor=gaussian):
    u, v = zeros_like(x), zeros_like(y)
    for p, (x_p, y_p) in enumerate(zip(x, y)):
        r2 = (x - x_p)**2 + (y - y_p)**2
        kf = bs_kernel_factor(r2, blob_size**2)
        K1, K2 = -kf * y, kf * x
        K1[p], K2[p] = 0.0, 0.0
        u[p], v[p] = dot(circ, K1), dot(circ, K2)
    return u, v

