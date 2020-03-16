"""
Combine data from multiple datasets with the programmable filter. You need
to click on both of them in the pipeline and then go 'Programmable Filter'

This calculates some mixing stuff, looking at vorticiy and things like that

@author: Nick
"""
print "Beginning imports"
from paraview.numpy_support import vtk_to_numpy
import vtk.numpy_interface.algorithms as algs
from numpy import zeros,eye,maximum

def get_array(name):
    for inp in inputs:
        if name in inp.CellData.keys():
            arr = inp.CellData[name]
            return arr
    else:
        raise KeyError('{} not found'.format(name))

def post_array(arr,name):
    print("Saving", name, "...")
    output.CellData.append(arr, name)
    return

print("Calculating vorticity")
u =  get_array('u')
v =  get_array('v')
w =  get_array('w')

nx = u.size
neq = 3

du = zeros((neq,neq,nx)) # du[i,j] = du_i / dx_j
du[0] = algs.gradient(u).T # Paraview's damn convention [nx,neq]
du[1] = algs.gradient(v).T
du[2] = algs.gradient(w).T

vort = zeros((neq,nx))

vort[0] = du[2,1] - du[1,2]
vort[1] = du[0,2] - du[2,0]
vort[2] = du[1,0] - du[0,1]
post_array(vort[0],'vort_x')
post_array(vort[1],'vort_y')
post_array(vort[2],'vort_z')

