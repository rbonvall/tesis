#!/usr/bin/env python

from __future__ import division
from numpy import *
from numpy.linalg import norm

def biot_savart_kernel(x, y):
    denominator = 2 * pi * (x**2 + y**2)
    return (-y/denominator, x/denominator)

def eval_velocity(x, y, circ):
    u, v = zeros_like(x), zeros_like(y)
    for p, (x_p, y_p) in enumerate(zip(x, y)):
        K1, K2 = biot_savart_kernel(x - x_p, y - y_p)
        K1[p], K2[p] = 0.0
        u[p], v[p] = dot(circ, K1), dot(circ, K2)
    return u, v


