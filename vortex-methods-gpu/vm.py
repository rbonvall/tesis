#!/usr/bin/env python

from __future__ import division
from numpy import *

_m6_funclist = [
    (-1/88)  * poly1d((1, -1)) * poly1d((60, -87, -87, 88, 88)),
    (1/176)  * poly1d((1, -1)) * poly1d((1, -2)) * poly1d((60, -261, 257, 68)),
    (-3/176) * poly1d((1, -2)) * poly1d((4, -17, 12)) * poly1d((1, -3))**2,
    poly1d((0,)),
]
def m6(x):
    x = absolute(x)
    condlist = [
        x < 1,
        logical_and(1 <= x, x < 2),
        logical_and(2 <= x, x < 3),
        x >= 3,
    ]
    return piecewise(x, condlist, _m6_funclist)


# interpolation kernel
def W(x, y):
    return m6(x) * m6(y)

def vorticity(u, v):
    ux, uy = gradient(u)
    vx, vy = gradient(v)
    return uy - vx


def remesh(w_p, x_p, y_p, x_mesh, y_mesh, h):
    mesh_width, mesh_height = x_mesh.shape
    #
    tiled_x_p = tile(x_p, (mesh_width, mesh_height, 1))
    tiled_y_p = tile(y_p, (mesh_width, mesh_height, 1))
    #
    reshaped_x_mesh = x_mesh.reshape(x_mesh.shape + (1,))
    reshaped_y_mesh = y_mesh.reshape(y_mesh.shape + (1,))
    #
    dx = reshaped_x_mesh - tiled_x_p
    dy = reshaped_y_mesh - tiled_y_p
    #
    h_inv = 1/h
    weights = W(h_inv * dx, h_inv * dy)
    w_mesh = dot(weights, w_p)
    return w_mesh


