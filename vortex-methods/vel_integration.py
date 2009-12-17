#!/usr/bin/env python

from numpy import zeros_like
import functools
import vm


def euler(x, y, circ, dt, squared_blob_size):
    u, v = vm.eval_velocity(x, y, circ=circ,
                            squared_blob_size=squared_blob_size)
    return u, v


def runge_kutta(x, y, circ, dt, squared_blob_size):
    u, v = zeros_like(x), zeros_like(y)

    eval_velocity = functools.partial(vm.eval_velocity, circ=circ,
                                      squared_blob_size=squared_blob_size)

    kx, ky = eval_velocity(x, y) # k1
    u += kx/6
    v += ky/6

    dx, dy = kx * (dt/2), ky * (dt/2)
    kx, ky = eval_velocity(x + dx, y + dy) # k2
    u += kx/3
    v += ky/3

    dx, dy = kx * (dt/2), ky * (dt/2)
    kx, ky = eval_velocity(x + dx, y + dy) # k3
    u += kx/3
    v += ky/3

    dx, dy = kx * dt, ky * dt
    kx, ky = eval_velocity(x + dx, y + dy) # k4
    u += kx/6
    v += ky/6

    return u, v

