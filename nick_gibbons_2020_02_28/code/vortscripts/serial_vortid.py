"""
Vortex visualisation using the Triple Decomposition Method

@author: Nick Gibbons
"""
from solvers import EEStwo
from paraview.numpy_support import vtk_to_numpy
import vtk.numpy_interface.dataset_adapter as dsa
import vtk.numpy_interface.algorithms as algs
from numpy import zeros,eye,minimum,sin,cos,radians,linspace,outer,save,convolve
from numpy.linalg import eigvals
from scipy.optimize import minimize, fmin
from time import time 


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

def ComputeQ(a,b,c):
    """ Compute the rotations matrix from kolar_vort07 appendix A """
    Q = zeros((3,3))
    Q[0,0] = cos(a)*cos(b)*cos(c) - sin(a)*sin(c)
    Q[0,1] = sin(a)*cos(b)*cos(c) + cos(a)*sin(c)
    Q[0,2] =                      - sin(b)*cos(c)
    Q[1,0] =-cos(a)*cos(b)*sin(c) - sin(a)*cos(c)
    Q[1,1] =-sin(a)*cos(b)*sin(c) + cos(a)*cos(c)
    Q[1,2] =                        sin(b)*sin(c)
    Q[2,0] = cos(a)*sin(b)
    Q[2,1] = sin(a)*sin(b)
    Q[2,2] = cos(b)
    return Q

def refscalar(angles, W, S):
    """ Compute the Reference Frame Scalar """
    a,b,c = angles
    Q = ComputeQ(a,b,c)
    Sdash = Q.dot(S)
    Sdash = Sdash.dot(Q.T)
    Wdash = Q.dot(W)
    Wdash = Wdash.dot(Q.T)
    REF = abs(Sdash[0,1]*Wdash[0,1])
    REF+= abs(Sdash[1,2]*Wdash[1,2])
    REF+= abs(Sdash[2,0]*Wdash[2,0])
    return -1*REF/1e6 # Negative so that we can use fmin

u = get_array('u')
v = get_array('v')
w = get_array('w')

nx = u.size
neq = 3
print "Starting with nx", nx

du = zeros((neq,neq,nx)) # du[i,j] = du_i / dx_j
du[0] = algs.gradient(u).T   # Paraview's damn convention [nx,neq]
du[1] = algs.gradient(v).T
du[2] = algs.gradient(w).T

W = zeros((neq,neq,nx))
S = zeros((neq,neq,nx))
L2= zeros(nx)

for i in range(neq):
    for j in range(neq):
        S[i,j] = 0.5*(du[i,j] + du[j,i])
        W[i,j] = 0.5*(du[i,j] - du[j,i])

print "du: "
print du[:,:,0]


t0 = time()
for n in range(nx):
    if n%100==0:
        print '@', round(float(n)/nx*100, 5), '%'
    Sn = S[:,:,n]
    Wn = W[:,:,n]

    print "S:"
    print Sn
    print "W:"
    print Wn
    guess = array([radians(20.0),radians(10.0), radians(-40.0)])
    result = refscalar(guess, Wn, Sn)
    print "Test result", result
    guess = array([radians(20.0),radians(10.0), radians(-40.0)])
    Q = ComputeQ(*guess)
    print "test Q", Q
    temp = Q.dot(Wn)
    Wdash = temp.dot(Q.T)
    print "test Wdash", Wdash


    result=fmin(refscalar, guess, (Wn, Sn), full_output=True,disp=False, xol=1e-6)
    angles, fopt, iter, fcalls, flag = result
    if flag!=0:
        print "Screwup with result, n: ", result, n
        raise Exception("Solver not converged")
    a,b,c = angles
    print "a,b,c", a,b,c
    # Now compute what?
    Q = ComputeQ(a,b,c)
    print "Q", Q
    dudash = (Q.dot(du[:,:,n])).dot(Q.T)
    resdash= dudash.copy()

    for i in range(3):
        for j in range(3):
            if i==j: continue
            resdash[i,j] = sign(dudash[i,j])*minimum(abs(dudash[i,j]), abs(dudash[j,i]))
    
    print "resdash", resdash
    sheardash = dudash - resdash
    
    res = ((Q.T).dot(resdash)).dot(Q)
    shear= ((Q.T).dot(sheardash)).dot(Q)
    
    el = zeros((3,3))
    rr = zeros((3,3))
    
    for i in range(3):
        for j in range(3):
            el[i,j] = 0.5*(res[i,j] + res[j,i])
            rr[i,j] = 0.5*(res[i,j] - res[j,i])
    
    LL = el.dot(el) + rr.dot(rr)
    vals = eigvals(LL)
    print "vals", vals
    l1, l2, l3 = sorted(vals)
    L2[n] = l2
    break
    if n==100000: break

t1 = time()
print "Time taken was: ", t1-t0," sec"
print "100000/nx", 100000/nx

post_array(du[2,1] - du[1,2],'vort_x')
post_array(du[0,2] - du[2,0],'vort_y')
post_array(du[1,0] - du[0,1],'vort_z')
post_array(L2,'L2')

