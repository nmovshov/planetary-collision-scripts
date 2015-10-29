#!/soft/scipy_0.13.0/CentOS_6/bin/python
#-------------------------------------------------------------------------------
# Quick and dirty plot of pressure vs. radius of nodes read from .fnl file(s).
#-------------------------------------------------------------------------------
import sys, os
import ahelpers

if len(sys.argv)==1:
    sys.exit("ERROR: provide file/directory name as first parameter.")

fdname = sys.argv[1]
if os.path.isfile(fdname):
    print "Plotting from file", fdname
    fnl = ahelpers.load_fnl(fdname)
    ahelpers.plot_P_vs_r(fnl, True)
elif os.path.isdir(fdname):
    print "Plotting from directory", fdname
    ahelpers.plot_P_vs_r_output(fdname, True)
else:
    sys.exit("ERROR: {} is not a valid file or directory.".format(fdname))
    pass
pass
