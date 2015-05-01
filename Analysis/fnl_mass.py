#!/soft/scipy_0.13.0/CentOS_6/bin/python
#---------------------------------------------------------------------------------
# fnl_mass - command-line utility the prints mass found in SPH output data. 
#
# Author: Naor Movshovitz (nmovshov at gee mail dot com)
#---------------------------------------------------------------------------------
import sys
import numpy as np
import argparse
#import ahelpers

cout = sys.stdout.write
parser = argparse.ArgumentParser()
parser.add_argument('filename', help="name of file containing node list data")
args = parser.parse_args()

# Load node list data
cout("Reading file...")
try:
    fnl = ahelpers.load_fnl(args.filename)
    cout("Done.\n")
    pos = np.vstack((fnl.x, fnl.y, fnl.z)).T
    vel = np.vstack((fnl.vx, fnl.vy, fnl.vz)).T
    m   = fnl.m
    print "Found {} nodes in {} node lists totaling {} kg.".format(
        fnl.nbNodes, np.unique(fnl.id).size, sum(m))
    for n in np.unique(fnl.id):
        print "    List {:g}: {:.6g} kg in {} nodes ({:.4g} kg/node).".format(
            n, sum(m[fnl.id == n]), sum(fnl.id == n),
            sum(m[fnl.id == n])/sum(fnl.id == n))
except StandardError:
    try:
        # FNL structure as it is in 2015-05-01
        raw = np.loadtxt(args.filename)
        cout("Done.\n")
        nl_id = raw[:,0]
        pos = raw[:,2:5]
        vel = raw[:,5:8]
        m   = raw[:,8]
        print "Found {} nodes in {} node lists totaling {} kg.".format(
            nl_id.size, np.unique(nl_id).size, sum(m))
        for n in np.unique(nl_id):
            print "    List {:g}: {:.6g} kg in {} nodes ({:.4g} kg/node).".format(
                n, sum(m[nl_id == n]), sum(nl_id == n),
                sum(m[nl_id == n])/sum(nl_id == n))
    except:
        raise StandardError("Could not read data from {}".format(
            args.filename))
