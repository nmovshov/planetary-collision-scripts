#!/soft/scipy_0.13.0/CentOS_6/bin/python
#-------------------------------------------------------------------------------
# Quick and dirty scatter plot of nodes read from .fnl file, in XY plane.
#-------------------------------------------------------------------------------
import sys, os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

if len(sys.argv)==1:
    sys.exit("ERROR: provide file name as first parameter.")

nodes = np.loadtxt(sys.argv[1])
if (nodes.ndim != 2) or (nodes.shape[1] != 15):
    sys.exit("{} does not appear to contain a valid flattened node list".format(
             sys.argv[1]))

print "Plotting nodes from file", sys.argv[1]
x = nodes[:,2]/1e3
y = nodes[:,3]/1e3
m = nodes[:,8]
rho = nodes[:,9]
h = nodes[:,13]/1e3

plt.figure()
plt.scatter(x, y, s=50, c=rho, alpha=0.8)
plt.xlabel('X [km]')
plt.ylabel('Y [km]')
plt.title(sys.argv[1])
plt.show()
print "Done."
