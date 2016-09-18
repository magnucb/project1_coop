import numpy as np
import astropy as ap
import matplotlib.pylab as pl
import sys
import os
import scipy.linalg as linalg
import time

"""
Solve (d/dx)^2 u = f for known f.
Double differentials using Euler forward and backward methods:
->(d/dx)^2 u = (u[i+1] - 2 u[i] + u[i-1])/h^2 = f
-> u[i+1] - 2 u[i] + u[i-1] = h^2*f = y[i]
-> generalize -> a[i]u[i-1] + b[i]u[i] + c[i]u[i+1] = y[i]
"""
#NB program defaults version to 0, if succesful run as diff. version.
curdir = os.getcwd()
try:
    n = int(sys.argv[1])
except IndexError:
    n = 5
except ValueError:
    sys.exit("Commandline-argument must be integer \nFirst number of gridpoints, then version number")

try:
    version = str(int(sys.argv[2]))
except IndexError:
    version = "0"
except ValueError:
    sys.exit("Commandline-argument must be integer \nFirst number of gridpoints, then version number")
    
def write2file(outstring,
               filename=curdir+"/data/testfile_v%d.dat",
               append=True):
    """
    If 'append' is True:
    -open 'filename'.
    -append 'outstring' to end of file
    -close file
    If 'append' is False:
    -open a new file 'filename' (deleting the old one)
    -write 'outstring' to file
    -close file
    """
    outstring = str(outstring)
    if append:
        with open(filename,'a') as outfile:
            outfile.write(outstring + "\n")
    else:
        with open(filename,'w') as outfile:
            outfile.write(outstring + "\n")
    return outstring

def general_tridiag(a, b, c, y):
    #args: arrays for tridiagonal matrix (below, on and above), array for vertical solution
    n = len(y) #number of matrix operations
    u = np.zeros(n) #first and last value must remain zero, dirichlet BC
    #forward substitution
    for i in range(1,n-1): #do not change vert[0] and vert[1]
        k = a[i]/float(b[i-1])
        b[i] -= k*c[i-1]
        y[i] -= k*y[i-1]
    #backward subtitution
    u[n-2] = y[n-2]/b[n-2] #second to last array-element
    for i in reversed(xrange(1,n-2)):
        u[i] = (y[i] - u[i+1]*c[i])/float(b[i])
    return u

def specific_tridiag(y):
    #args: arrays for tridiagonal matrix (below, on and above), array for vertical solution
    n = len(y) #number of matrix operations
    u = np.zeros(n) #first and last value must remain zero, dirichlet BC
    d = np.array([(i+1)/float(i) for i in range(1,n)])
    #forward substitution
    for i in range(1,n-1): #do not change vert[0] and vert[1]
        y[i] = y[i] - y[i-1]/float(d[i-1])
    #backward subtitution
    u[n-2] = y[n-2]/d[n-2] #second to last array-element
    for i in reversed(xrange(1,n-2)):
        u[i] = (y[i] + u[i+1])/float(d[i])
    return u

def general_LU_decomp(A_matrix, vert):
    if len(vert)>1000:
        return np.array([False])
    #arg: matrix of linear equation, array of solution
    u = np.zeros(len(vert))
    A_matrix = A_matrix[1:-1,1:-1]
    vert = vert[1:-1] #do not work on Dirichlet BC
    P, L, U = linalg.lu(a = A_matrix, overwrite_a = False, check_finite = True) #scipy-function
    # solve L*w = vert. #overwrite the array vert.
    w = linalg.solve(U,vert, sym_pos=True, lower=True)
    # solve L*w = vert. #overwrite the array vert.
    u[1:-1] = linalg.solve(U,vert, sym_pos=True, lower=True)
    return u

def test_diag(d):
    #d must be 1-D array
    n = len(d)
    exact = np.zeros(n)
    for i in range(n):
        exact[i] = - float(i+2)/float(i+1)
    return abs((d-exact)/exact)

x0=0.0; x1=1.0; h=(x1-x0)/(n+1.0)
#vectors of tridiagonal matrix A
a = -1.0*np.ones(n)
b = 2.0*np.ones(n)
c = -1.0*np.ones(n)
A = np.diag(a[1:], -1) + np.diag(b, 0) + np.diag(c[:-1], 1)
#vectors y and x
x = np.linspace(x0,x1,n) 
y = h*h*100*np.exp(-10*x)
u_exact = 1-(1-np.exp(-10))*x-np.exp(-10*x)

t0 = time.clock()

#general gaussian elimination of tri-diagonal matrix
u_gen = general_tridiag(a, b, c, y)
t1 = time.clock()
t_gen = t1 - t0

#specific gaussian elimination of tri-diagonal matrix
u_spec = specific_tridiag(y)
t2 = time.clock()
t_spec = t2 - t1

#LU-decomposition of tri-diagonal matrix
u_LU = general_LU_decomp(A_matrix=A, vert=y)
t3 = time.clock()
t_LU = t3 - t2

"""
#double check diagonal of matrix
diag_error = test_diag(b)
max_diag_error = max(diag_error)
if max_test_diag >= 1e-14:
    print "maximum error in diagonal is: e=", max_test_diag
"""

#storing data
datafile_u = curdir+"/data/dderiv_u_python_v%s_n%d.dat"%(version,n) #datafile for calculated arrays
write2file("calculated arrays for u(x) for n=%d"%n, filename=datafile_u, append=False) # create first line of datafile
datafile_time = curdir+"/data/dderiv_time_python_v%s_n%d.dat"%(version,n) #datafile for calculated arrays
write2file("CPU time for the diffenrent", filename=datafile_time, append=False) #create first line of datafile

if not u_LU.any():
    #store data for u(x) for three different methods
    write2file("x, f, u_gen, u_spec", filename=datafile_u, append=False) # create first line of datafile
    for i in range(len(x)):
        write2file("%1.10e, %1.10e, %1.10e, %1.10e"%(x[i], -y[i]/h/h, u_gen[i], u_spec[i]),
                   filename=datafile_u,
                   append=True)
    #store data for CPU time
    datastring_base = "CPU-time for method: "
    methods = ["general_tridiagonal", "specific_tridiagonal"]
    CPUtime = [t_gen, t_spec]
    for i in range(len(methods)):
        print write2file(datastring_base+"%s, %1.10e s"%(methods[i], CPUtime[i]),
                         filename=datafile_time,
                         append=True)

else:
    #store data for u(x) for three different methods
    write2file("x, f, u_gen, u_spec, u_LU", filename=datafile_u, append=False) # create first line of datafile
    for i in range(len(x)):
        write2file("%1.10e, %1.10e, %1.10e, %1.10e, %1.10e"%(x[i], -y[i]/h/h, u_gen[i], u_spec[i], u_LU[i]),
                   filename=datafile_u,
                   append=True)
    #store data for CPU time
    datastring_base = "CPU-time for method: "
    methods = ["general_tridiagonal", "specific_tridiagonal", "LU-decomposition"]
    CPUtime = [t_gen, t_spec, t_LU]
    for i in range(len(methods)):
        print write2file(datastring_base+"%s, %1.10e s"%(methods[i], CPUtime[i]),
                         filename=datafile_time,
                         append=True)
