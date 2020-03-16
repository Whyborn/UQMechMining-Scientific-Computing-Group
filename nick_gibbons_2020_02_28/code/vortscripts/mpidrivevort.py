"""
This is a sort of control program to drive the TDM algorithm in MPI mode,
using python to iterface with paraview in a relatively simple way.

@author: Nick Gibbons
"""

from __future__ import division, print_function
from mpi4py import MPI # calls MPI_Init()
from numpy import save, load, array, zeros, concatenate
from ctypes import cdll,c_double,POINTER,c_int
from time import time

mcw = MPI.COMM_WORLD
rank = mcw.Get_rank()
ws = mcw.Get_size()

for i in range(ws):
    if i==rank:
        print("Rank i loading", i)
        du = load('du.npy')
        nx,neq,neq = du.shape
        
        # Load balance problem AKA ceil division
        if nx%ws==0: chunksize = nx//ws
        else: chunksize = nx//ws + 1
        
        srt = chunksize*rank
        end = chunksize*(rank+1)
        mydu= du[srt:end].copy()
        
        del du
    mcw.Barrier()

mynx, myneq, myneq = mydu.shape
print(srt, end, mydu.shape, mynx)
myL = zeros((mynx))
if rank==0: verbose=1
else: verbose=0

# Launch the core code from libvortid.so
lib = cdll.LoadLibrary('../scripts/vortscripts/libvortid.so')
c_double_p = POINTER(c_double)
du_c = mydu.ctypes.data_as(c_double_p)
L_c  = myL.ctypes.data_as(c_double_p)

# We're allowed to pass ints
print("Launching Compute Core on Rank {}".format(rank))
lib.ComputeL2TDM(du_c, mynx, L_c, verbose) 
#lib.ComputeL2(du_c, mynx, L_c, verbose) 
print("Rank {} Done Compute...".format(rank))

# Sync data and rebuild the array on rank 0
if rank==0:
    L2s = [myL]
    for i in range(1,ws):
        othernxbuffer = array([0])
        mcw.Recv([othernxbuffer, MPI.INT], source = i, tag=i)
        othernx = othernxbuffer[0]
        print("Got othernx of {}".format(othernx))

        otherbuffer = zeros((othernxbuffer[0]))
        mcw.Recv([otherbuffer, MPI.DOUBLE], source = i, tag=i)
        L2s.append(otherbuffer)

    L2 = concatenate(L2s)
    print("Saving {} length L2...".format(L2.size))
    save('L2TDM_mean.npy', L2)

else:
    size = array([mynx])
    mcw.Send([size,MPI.INT], dest=0, tag=rank)
    mcw.Send([myL,MPI.DOUBLE], dest=0, tag=rank)

print("Rank {} Done. Call Shutdown".format(rank))
