#! /proj/nmovshov_hindmost/collisions/SPHERAL/bin/python
#-------------------------------------------------------------------------------
# Quick and dirty plot of pressure vs. radius of nodes read from .fnl file.
#-------------------------------------------------------------------------------
import sys
import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import shelpers

if len(sys.argv)==1:
    sys.exit("Please provide file name as first parameter")

nodes = sp.loadtxt(sys.argv[1])
if (nodes.ndim != 2) or (nodes.shape[1] != shelpers.nb_fnl_columns):
    sys.exit("{} does not appear to contain a valid flattened node list".format(
             sys.argv[1]))

print "Plotting nodes from file", sys.argv[1]
x = nodes[:,2]
y = nodes[:,3]
z = nodes[:,4]
r = np.sqrt(x**2 + y**2 + z**2)
P = nodes[:,10] * 1e-9 # pressure in GPa

plt.figure()
plt.plot(r,P,'.')
plt.xlabel('Radius [m]')
plt.ylabel('Pressure [GPa]')
plt.title(sys.argv[1])
plt.grid()
plt.show()
print "Done."
