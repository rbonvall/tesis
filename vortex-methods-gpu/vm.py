#!/usr/bin/env python

from __future__ import division
from numpy import *
import scipy.signal
import scipy.sparse
import scipy.linalg.iterative
from math import ceil

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
    uy, ux = gradient(u)
    vy, vx = gradient(v)
    return vx - uy

def poisson_source(w_m):
    w_y, w_x = gradient(w_m)
    return -w_y, w_x  # = -curl(w)

def poisson_solve(f, g, h):
    M, N = f.shape
    # sparse CSC representation of FDM matrix for Poisson eq.
    d0 = 4 * ones(M * N)
    d1 =    -ones(M * N); d1[M-1::M] = 0
    d2 =    -ones(M * N)
    A = scipy.sparse.spdiags([d2, d1, d0, d1, d2],
                             [-M, -1,  0,  1,  M], M * N, M * N)

    # Solve first component of equation (lapl(u) = -w_y)

    # boundaries
    # (assumption: periodic boundary conditions)
    # TODO: handle other kinds of BCs
    top_boundary    = f[:,  0]
    bottom_boundary = f[:, -1]
    left_boundary   = f[-1, :]
    right_boundary  = f[ 0, :]

    # build right-hand side vector of the FDM system
    b = f.reshape(f.size).copy() * h**2
    #b[       :M] -= bottom_boundary
    #b[(N-1)*M: ] -= top_boundary
    #b[   ::M] -= left_boundary
    #b[M-1::M] -= right_boundary

    # solve for each velocity component using conjugate gradient method
    u0 = b/4
    u, u_info = scipy.linalg.iterative.cg(A, b, x0=u0, tol=1e-8)
    if u_info != 0:
        raise RuntimeError("Poisson solver did not converge")

    subplot(131)
    plot(b)
    subplot(132)
    plot(u0)
    subplot(133)
    plot(u)


    # Solve second component of equation (lapl(v) = w_x)

    # boundaries
    # (assumption: periodic boundary conditions)
    # TODO: handle other kinds of BCs
    top_boundary    = g[:,  0]
    bottom_boundary = g[:, -1]
    left_boundary   = g[-1, :]
    right_boundary  = g[ 0, :]

    # build right-hand side vector of the FDM system
    b = g.reshape(g.size).copy() * h**2
    #b[       :M] -= bottom_boundary
    #b[(N-1)*M: ] -= top_boundary
    #b[   ::M] -= left_boundary
    #b[M-1::M] -= right_boundary

    # solve for each velocity component using conjugate gradient method
    v0 = b/4
    v, v_info = scipy.linalg.iterative.cg(A, b, x0=v0, tol=1e-8)
    if v_info != 0:
        raise RuntimeError("Poisson solver did not converge")


    return u.reshape(f.shape), v.reshape(g.shape)


def BROKEN_remesh(w_p, x_p, y_p, x_mesh, y_mesh, h):
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


def SLOW_remesh(w_p, x_p, y_p, x_m, y_m, h):
    #mesh_width, mesh_height = x_m.shape
    nr_particles = w_p.size

    w_m = zeros_like(x_m)
    h_inv = 1/h
    for p in xrange(nr_particles):
        if p % 256 == 0: print "   ", p
        dx = x_m - x_p[p]
        dy = y_m - y_p[p]
        weights = W(dx * h_inv, dy * h_inv)
        w_m += w_p[p] * weights
    return w_m

    
def remesh(w_p, x_p, y_p, x_m, y_m, h):
    M, N = x_m.shape
    nr_particles = w_p.size

    w_m = zeros_like(x_m)
    h_inv = 1/h
    for p in xrange(nr_particles):
        print p

        # TODO: dub magic constants in terms of interpolation kernel's support
        i_p, j_p = x_p[p]/h, y_p[p]/h
        i_c, j_c = ceil(i_p), ceil(j_p)
        di, dj = i_c - i_p, j_c - j_p
        i_int, j_int = mgrid[di-3:di+3, dj-3:dj+3]
        weights = W(i_int, j_int)[max(0, 3 - i_c):6 - max(0, i_c + 3 - M),
                                  max(0, 3 - j_c):6 - max(0, j_c + 3 - N)]
        w_m[max(0, i_c - 3):min(i_c + 3, M),
            max(0, j_c - 3):min(j_c + 3, N)] += w_p[p] * weights
    return w_m



def diffusion(w_m, h):
    diffusion_stencil = array([[0,  1, 0],
                               [1, -4, 1],
                               [0,  1, 0]], dtype=float) / h**2
    return scipy.signal.correlate2d(w_m, diffusion_stencil,
                                    mode='same', boundary='wrap')



