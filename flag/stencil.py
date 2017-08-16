#!/usr/bin/env python3
# ----------------------------------------- #
#
import numpy as np
from scipy.linalg import solve

def get_ghost(h,known_vals,known_powers,unknown_powers,y):
    y0 = np.inner(known_vals,[h**i for i in known_powers])
    A = np.matrix([[h**i for i in unknown_powers] for h in [0.5,1.5,2.5]])
    b = np.array(y[0:3])-y0
    unknown_vals = solve(A,b)
    return y0+np.inner(unknown_vals,[h**i for i in unknown_powers])

def pinned_ghost(h,y0,y):
    return get_ghost(h,[y0,0],[0,2],[1,3,4],y)

def clamped_ghost(h,y0,y):
    return get_ghost(h,[y0,0],[0,1],[2,3,4],y)

def free_ghost(h,y0,y):
    return get_ghost(h,[0,0],[2,3],[0,1,4],y)

def get_4th(i,y0,y,ghost_fun):
    y_ghost = np.concatenate(([ghost_fun(h,y0,y) for h in [-1.5,-0.5]],y))
    j = i+2
    return y_ghost[j+2]-4*y_ghost[j+1]+6*y_ghost[j]-4*y_ghost[j-1]+y_ghost[j-2]

def get_stencil(i,ghost_fun):
    A = np.identity(6)
    return [get_4th(i,A[j,0],A[j,1:],ghost_fun) for j in range(6)]

import unittest
class TestMethods(unittest.TestCase):
    def test_ghost(self):
        self.assertAlmostEqual(pinned_ghost(-0.5,1,[1.5,2.5,3.5]),0.5)
        self.assertAlmostEqual(free_ghost(-0.5,0,[0.5,1.5,2.5]),-0.5)
        self.assertAlmostEqual(clamped_ghost(-0.5,0,[0.5**2,1.5**2,2.5**2]),0.25)

    def test_4th(self):
        y0 = 2
        y = [y0+(x+0.5)**3 for x in range(4)]
        self.assertAlmostEqual(get_4th(0,y0,y,pinned_ghost),0)
        y = [y0+(x+0.5)**4/24 for x in range(4)]
        self.assertAlmostEqual(get_4th(1,y0,y,pinned_ghost),1)
        y = [y0+(x+0.5) for x in range(4)]
        self.assertAlmostEqual(get_4th(0,y0,y,free_ghost),0)
        y = [y0+(x+0.5)**4/24 for x in range(4)]
        self.assertAlmostEqual(get_4th(1,y0,y,free_ghost),1)
        y = [y0+(x+0.5)**3 for x in range(4)]
        self.assertAlmostEqual(get_4th(0,y0,y,clamped_ghost),0)
        y = [y0+(x+0.5)**4/24 for x in range(4)]
        self.assertAlmostEqual(get_4th(1,y0,y,clamped_ghost),1)

    def test_stencil(self):
        self.assertAlmostEqual(get_stencil(2,clamped_ghost),[0,1,-4,6,-4,1])

# unittest.main()

print('pinned')
print(get_stencil(0,pinned_ghost))
print(get_stencil(1,pinned_ghost))

print('clamped')
print(get_stencil(0,clamped_ghost))
print(get_stencil(1,clamped_ghost))

print('free')
print(get_stencil(0,free_ghost))
print(get_stencil(1,free_ghost))
