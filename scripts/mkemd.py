#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import re
import sys
import numpy as np
import itertools

class FileParser:
    def __init__(self, fn):
        s = ""
        f = open(fn)

        with open(fn) as f:
            for line in f:
                s += re.sub(r"#.*$", "", line)
        self.numbers = s.split()

    def next_int(self):
        return int(self.numbers.pop(0))

    def next_float(self):
        return float(self.numbers.pop(0))

class EmdGenerator:
    def read_pmt(self, fn):
        fp = FileParser(fn)

        nfreq = fp.next_int()
        self.freq = []
        for i in range(0, nfreq):
            self.freq.append(fp.next_float())

        nobsgrp = fp.next_int()
        obsgrps = []
        for i in range(0, nobsgrp):
            nobstypeidx = fp.next_int()
            obstypeidxes = np.zeros(nobstypeidx, dtype = 'int32')
            for j in range(0, nobstypeidx):
                obstypeidxes[j] = fp.next_int()

            nfreqidx = fp.next_int()
            freqidxes = np.zeros(nfreqidx, dtype = 'int32')
            for j in range(0, nfreqidx):
                freqidxes[j] = fp.next_int()

            rxrange = np.zeros((2, 3), dtype = 'float32')
            for j in range(0, 2):
                for k in range(0, 3):
                    rxrange[j][k] = fp.next_float()

            obsgrps.append((obstypeidxes, freqidxes, rxrange))

        rxes = set()
        for i in range(0, len(obsgrps)):
            rxrange = obsgrps[i][2]
            y = np.arange(rxrange[0][0], rxrange[0][1] + 0.001, rxrange[0][2])
            z = np.arange(rxrange[1][0], rxrange[1][1] + 0.001, rxrange[1][2])
            for rx in itertools.product(*[y, z]):
                rxes.add(rx)

        rxes = sorted(rxes)
        self.rxtoidx = {}
        self.idxtorx = {}
        idx = 0
        for rx in rxes:
            self.rxtoidx[rx] = idx
            self.idxtorx[idx] = rx
            idx += 1

        self.obses = set()
        for obsgrp in obsgrps:
            rxrange = obsgrp[2]
            y = np.arange(rxrange[0][0], rxrange[0][1] + 0.001, rxrange[0][2])
            z = np.arange(rxrange[1][0], rxrange[1][1] + 0.001, rxrange[1][2])
            rxes = itertools.product(*[y, z])
            rxidxes = []
            for rx in rxes:
                rxidxes.append(self.rxtoidx[rx])
            for obs in itertools.product(*[obsgrp[0], obsgrp[1], rxidxes]):
                self.obses.add(obs)

        self.obses = sorted(self.obses)

    def save_emd(self, fn):
        emdf = open(fn, 'w')

        print("# frequencies", file = emdf)
        print("%d" % (len(self.freq)), file = emdf)
        for f in self.freq:
            print("% .6E" % (f), file = emdf)

        print("# receivers", file = emdf)
        print("%d" % (len(self.idxtorx)), file = emdf)
        for i in range(0, len(self.idxtorx)):
            rx = self.idxtorx[i]
            print("% .6E % .6E" % (rx[0], rx[1]), file = emdf)

        print("# observations", file = emdf)
        print("%d" % (len(self.obses)), file = emdf)
        for o in self.obses:
            print("%7d %7d %7d % .6E % .6E % .6E % .6E" % (o[0], o[1], o[2], 1.0, 1.0, 1.0, 1.0), file = emdf)

if __name__ == '__main__' :
    if len(sys.argv) < 2 :
        print("Parameter is NOT enough.")
        sys.exit()

    eg = EmdGenerator()
    eg.read_pmt(sys.argv[1] + '-emd.pmt')
    eg.save_emd(sys.argv[1] + '.emd')
