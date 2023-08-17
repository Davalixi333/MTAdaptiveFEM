#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import print_function

import re
import sys
import numpy as np

class FileParser:
    def __init__(self, fn):
        s = ""

        with open(fn) as f:
            for line in f:
                s += re.sub(r"#.*$", "", line)
        self.numbers = s.split()

    def next_int(self):
        return int(self.numbers.pop(0))

    def next_float(self):
        return float(self.numbers.pop(0))

class RSPConverter:
    def read_rsp(self, fn):
        fp = FileParser(fn)

        nfreq = fp.next_int()
        freqs = []
        for i in range(0, nfreq):
            freqs.append(fp.next_float())

        nrx = fp.next_int()
        rxes = np.zeros((nrx, 2), dtype = 'float32')
        for i in range(0, nrx):
            rxes[i][0] = fp.next_float()
            rxes[i][1] = fp.next_float()

        nobs = fp.next_int()
        self.obses = np.zeros((nobs, 10), dtype = 'float32')
        for i in range(0, nobs):
            self.obses[i][0] = fp.next_int()
            fidx = fp.next_int()
            tidx = fp.next_int()
            ridx = fp.next_int()
            fp.next_int()

            self.obses[i][1] = freqs[fidx]
            self.obses[i][2] = tidx
            self.obses[i][3] = rxes[ridx][0]
            self.obses[i][4] = rxes[ridx][1]
            self.obses[i][5] = fp.next_float()
            self.obses[i][6] = fp.next_float()
            self.obses[i][7] = fp.next_float()
            self.obses[i][8] = fp.next_float()

    def save_dat(self, fn):
        datf = open(fn, "w")

        obses = self.obses
        for i in range(obses.shape[0]):
            print("%2d %.6E % 2d % .5E % .5E % .6E % .6E % .6E % .6E"
                % (obses[i][0], obses[i][1], obses[i][2], obses[i][3], obses[i][4],
                obses[i][5], obses[i][6], obses[i][7], obses[i][8]), file=datf)

if __name__ == '__main__' :
    if len(sys.argv) < 2 :
        print("Parameter is NOT enough.")
        sys.exit()

    rc = RSPConverter()
    rc.read_rsp(sys.argv[1])
    rc.save_dat(sys.argv[1].replace('rsp', 'dat'))
