#!/usr/bin/env python
# vim: set fileencoding=utf-8

import optparse
import sys
from numpy import *

def main():
    parser = optparse.OptionParser()
    (options, args) = parser.parse_args()
    if len(args) < 2:
        sys.stderr.write('Usage: datafile1 datafile2\n')
        sys.exit(-1)
    f1, f2 = args[:2]

    data1 = loadtxt(f1)
    data2 = loadtxt(f2)
    u1, v1 = data1[:, 3], data1[:, 4]
    u2, v2 = data2[:, 3], data2[:, 4]

    du_sq = (u1 - u2) ** 2
    dv_sq = (v1 - v2) ** 2

    e_p = sqrt((u1 - u2)**2 + (v1 - v2)**2)
    error = sum(e_p)
    print  "# nr_particles: {0}".format(len(e_p))
    print u"# error: Σ √(Δu² + Δv²):"
    print error


if __name__ == '__main__':
    main()




