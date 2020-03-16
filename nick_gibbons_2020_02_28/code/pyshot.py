"""
Test stub for showing off c integration

To run:
 $ gcc -c -fPIC odeint.c
 $ gcc -shared odeint.o -o libodeint.so

@author: Nick Gibbons
"""

from pylab import plot,show,title,xlabel,ylabel
from numpy import array, zeros
import time

dt = 0.0001
N = 61500
x  = array([0.0, 0.0])
xd = array([30.0, 30.0])
xdd= array([0.0, -9.81])
xs = zeros((N,2))

start = time.clock()
for i in range(N):
    x += xd*dt + 0.5*xdd*dt**2
    xd += xdd*dt
    xs[i] = x.copy()

end = time.clock()
print("Time for N={}: {:4.4f}ms".format(N, (end-start)*1000.0))

plot(xs[:,0], xs[:,1])
title("Pure Python Particle Motion Simulation")
xlabel('X (m)')
ylabel('Y (m)')
show()

