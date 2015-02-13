#---------------------------------------------------------------------------------
# fnl_mass - command-line utility the prints mass found in SPH output data. 
#
# Author: Naor Movshovitz (nmovshov at gee mail dot com)
#---------------------------------------------------------------------------------
import sys
import numpy as np
import argparse
import ahelpers

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
        raw = np.loadtxt(args.filename)
        cout("Done.\n")
        pos = raw[:,0:3]
        vel = raw[:,3:6]
        m   = raw[:,6]
        print "Found {} particles totaling {} kg.".format(
            len(pos),sum(m))
    except:
        raise StandardError("Could not read data from {}".format(
            args.filename))
