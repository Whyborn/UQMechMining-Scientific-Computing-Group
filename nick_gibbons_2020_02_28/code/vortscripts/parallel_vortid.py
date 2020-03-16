"""
Vortex visualisation using the Triple Decomposition Method
This is the parallel version that uses some next level 
wizradry to save the data to disk and call a separate script
to process it.

@author: Nick Gibbons
"""
from paraview.numpy_support import vtk_to_numpy
import vtk.numpy_interface.dataset_adapter as dsa
import vtk.numpy_interface.algorithms as algs
from numpy import zeros,eye,minimum,sin,cos,radians,linspace,outer,save,convolve
from time import time 
from subprocess import call 


print("Begin...")
def get_array(name):
    for inp in inputs:
        if name in inp.CellData.keys():
            arr = inp.CellData[name]
            return arr
    else:
        raise KeyError('{} not found'.format(name))

def numpy2vtk(arr,dset,aa):
    vtkdata = dsa.numpyTovtkDataArray(arr)
    vtkarray = dsa.vtkDataArrayToVTKArray(vtkdata,dset)
    vtkarray.Association = aa
    return vtkarray

def post_array(arr,name):
    print("Saving", name, "...")
    output.CellData.append(arr, name)
    return

u = get_array('um')
v = get_array('vm')
w = get_array('wm')

nx = u.size
neq = 3
print "Starting with nx", nx

du = zeros((nx,neq,neq)) # du[0,i,j] = du_i / dx_j
du[:,0,:] = algs.gradient(u)  
du[:,1,:] = algs.gradient(v)
du[:,2,:] = algs.gradient(w)

save('du.npy', du)
del du

t0 = time()
call(['mpirun','-np','4','python3','../scripts/vortscripts/mpidrivevort.py'])
t1 = time()
print "Time taken was: ", t1-t0," sec"

#L2TDM = load('L2TDM.npy')
#L2 = load('L2.npy')

#post_array(du[:,2,1] - du[:,1,2],'vort_x')
#post_array(du[:,1,0] - du[:,0,1],'vort_z')
#post_array(L2TDM,'L2TDM')
#post_array(L2,'L2')

