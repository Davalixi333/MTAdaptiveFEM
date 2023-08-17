#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function

import re
import sys
import datetime
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

def gen_grid(grid, start):
    d = []
    coord = []

    for i in range(0, len(grid)):
        block = []
        for j in range(0, grid[i][2]):
            block.append(int(grid[i][0] * (grid[i][1] ** (j + 1))))
        if grid[i][3] == -1:
            d.extend(reversed(block))
        elif grid[i][3] == 1:
            d.extend(block)
        else:
            pass

    coord.append(start)
    for x in d:
        coord.append(coord[-1] + x)

    return (np.array(d), np.array(coord))

class ModGenerator:
    def read_pmt(self, fn):
        fp = FileParser(fn)

        ygrid = []
        ny = fp.next_int()
        ystart = fp.next_float()
        for i in range(0, ny):
            ygrid.append((fp.next_float(), fp.next_float(), fp.next_int(), fp.next_int()))

        zgrid = []
        nz = fp.next_int()
        zstart = fp.next_float()
        for i in range(0, nz):
            zgrid.append((fp.next_float(), fp.next_float(), fp.next_int(), fp.next_int()))

        blocks = []

        nlay = fp.next_int()
        for i in range(0, nlay):
            blocks.append((-1E15, 1E15, fp.next_float(), fp.next_float(), fp.next_int()))

        nblk = fp.next_int()
        for i in range(0, nblk):
            blocks.append((fp.next_float(), fp.next_float(), fp.next_float(), fp.next_float(), fp.next_int()))

        nrho = fp.next_int()

        self.rho = np.zeros(nrho, dtype='float64')
        self.rho_name = ['rho']

        self.lb = np.zeros(nrho, dtype='float64');
        self.ub = np.zeros(nrho, dtype='float64');

        for i in range(0, nrho):
            self.rho[i] = fp.next_float()
            self.lb[i] = fp.next_float()
            self.ub[i] = fp.next_float()

        self.dy, self.ycoord = gen_grid(ygrid, ystart)
        self.dz, self.zcoord = gen_grid(zgrid, zstart)

        self.cell_idx = np.zeros((len(self.dy), len(self.dz)), dtype = 'int32')

        self.insert_blk(blocks)

    def insert_blk(self, blocks):
        for b in range(0, len(blocks)) :
            for j in range(0, len(self.dy)) :
                y = self.ycoord[j] + self.dy[j] * 0.5
                for k in range(0, len(self.dz)) :
                    z = self.zcoord[k] + self.dz[k] * 0.5
                    if y > blocks[b][0] and y < blocks[b][1] and z > blocks[b][2] and z < blocks[b][3] :
                        self.cell_idx[j][k] = blocks[b][4]

    def save_mdl(self, fn):
        npoints = len(self.ycoord) * len(self.zcoord)
        points = np.zeros((npoints, 2), dtype='float64')

        p = 0
        for k in range(0, len(self.zcoord)):
            for j in range(0, len(self.ycoord)):
                    points[p][0] = self.ycoord[j]
                    points[p][1] = self.zcoord[k]
                    p = p + 1

        ncells = len(self.dy) * len(self.dz)
        cells = np.zeros((ncells, 4), dtype='int32')
        rho = np.zeros(ncells, dtype='float64')
        lb = np.zeros(ncells, dtype='float64')
        ub = np.zeros(ncells, dtype='float64')

        ny = len(self.dy)
        nz = len(self.dz)
        for z in range(0, nz):
            for y in range(0, ny):
                c = y + z * ny
                cells[c][0] = y + z * (ny + 1)
                cells[c][1] = y + z * (ny + 1) + 1
                cells[c][2] = y + (z + 1) * (ny + 1)
                cells[c][3] = y + (z + 1) * (ny + 1) + 1
                rho[c] = self.rho[self.cell_idx[y][z]]
                lb[c] = self.lb[self.cell_idx[y][z]]
                ub[c] = self.ub[self.cell_idx[y][z]]

        with open(fn + '.mdl', 'w') as mdlf:
            print('# ycoord dy', file=mdlf)
            print('# %d' % (ny), file=mdlf)
            for i in range(0, ny):
                print('# % E %E' % (self.ycoord[i], self.dy[i]), file=mdlf)
            print('# % E' % (self.ycoord[i + 1]), file=mdlf)
            print('', file=mdlf)

            print('# zcoord dz', file=mdlf)
            print('# %d' % (nz), file=mdlf)
            for i in range(0, nz):
                print('# % E %E' % (self.zcoord[i], self.dz[i]), file=mdlf)
            print('# % E' % (self.zcoord[i + 1]), file=mdlf)
            print('', file=mdlf)

            print('# points', file=mdlf)
            print('%d' % (npoints), file=mdlf)
            for p in range(0, npoints):
                print('% E % E' % (points[p][0], points[p][1]), file=mdlf)
            print('', file=mdlf)

            print('# cells', file=mdlf)
            print('%d 4' % (ncells), file=mdlf)
            for c in range(0, ncells):
                print('%8d %8d %8d %8d %8d' % (cells[c][0], cells[c][1], cells[c][2], cells[c][3], c), file=mdlf)
            print('', file=mdlf)

            print('# attributes', file=mdlf)
            print('%d' % (ncells), file=mdlf)
            for c in range(0, ncells):
                print('%E %E %E' % (rho[c], lb[c], ub[c]), file=mdlf)
            print('', file=mdlf)

        with open(fn + '.vtk', 'w') as vtkf:
            print('# vtk DataFile Version 3.0', file=vtkf)
            print('# This file was generated by mkmdl.py at ' + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), file=vtkf)
            print('ASCII', file=vtkf)
            print('DATASET UNSTRUCTURED_GRID', file=vtkf)
            print('', file=vtkf)

            print('POINTS %d double' % (npoints), file=vtkf)
            for p in range(0, npoints):
                print('% E % E 0' % (points[p][0], points[p][1]), file=vtkf)
            print('', file=vtkf)

            print('CELLS %d %d' % (ncells, ncells * 5), file=vtkf)
            for c in range(0, ncells):
                print('4 %8d %8d %8d %8d' % (cells[c][0], cells[c][1], cells[c][2], cells[c][3]), file=vtkf)
            print('', file=vtkf)

            print('CELL_TYPES %d' % (ncells), file=vtkf)
            for c in range(0, ncells):
                print('8', file=vtkf)
            print('', file=vtkf)

            print('CELL_DATA %d' % (ncells), file=vtkf)
            print('SCALARS %s double 1\nLOOKUP_TABLE default' % ('rho'), file=vtkf)
            for c in range(0, ncells):
                print('%E' %(rho[c]), file=vtkf)

if __name__ == '__main__' :
    if len(sys.argv) < 2 :
        print("Parameter not enough.")
        sys.exit()

    mg = ModGenerator()
    mg.read_pmt(sys.argv[1] + '-mdl.pmt')
    mg.save_mdl(sys.argv[1])
