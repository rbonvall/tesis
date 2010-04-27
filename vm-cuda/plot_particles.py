#!/usr/bin/env python

from __future__ import division
from numpy import *
import pylab
import sys

data = loadtxt(sys.stdin)
x = data[:, 0]
y = data[:, 1]
c = data[:, 2]
u = data[:, 3]
v = data[:, 4]

pylab.scatter(x, y, s=20, c=c, edgecolors='none')
pylab.quiver(x, y, u, v, width=0.001)
pylab.show()

