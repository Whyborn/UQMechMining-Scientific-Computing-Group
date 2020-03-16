"""
Test stub for showing off c integration

@author: Nick Gibbons
"""

from pylab import plot,show
from numpy import array, zeros
from ctypes import cdll,c_double,POINTER,c_int
import time

dt = 0.0001
N = 61500
x  = array([0.0, 0.0])
xd = array([30.0, 30.0])
xdd= array([0.0, -9.81])
xs = zeros((N,2))

c_double_p = POINTER(c_double)
xp = x.ctypes.data_as(c_double_p)
xdp = xd.ctypes.data_as(c_double_p)
xddp = xdd.ctypes.data_as(c_double_p)
xsp = xs.ctypes.data_as(c_double_p)
dt_c = c_double(dt)
lib = cdll.LoadLibrary('libodeint.so')

start = time.clock()
lib.integrate_particle(xp, xdp, xddp, xsp, N, 2, dt_c)
end = time.clock()

print("Time for N={}: {:4.4f}ms".format(N, (end-start)*1000.0))

plot(xs[:,0], xs[:,1])
show()

