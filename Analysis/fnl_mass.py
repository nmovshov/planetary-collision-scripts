#!/soft/scipy_0.13.0/CentOS_6/bin/python
#---------------------------------------------------------------------------------
# fnl_mass - command-line utility the prints mass found in SPH output data. 
#
# Author: Naor Movshovitz (nmovshov at gee mail dot com)
#---------------------------------------------------------------------------------
import sys
import numpy as np
import argparse

cout = sys.stdout.write
parser = argparse.ArgumentParser()
parser.add_argument('filename', help="name of file containing node list data")
parser.add_argument('-pr', '--print_radii', action='store_true',
    help="attempt to identify mean radii, assuming spherical geometry")
args = parser.parse_args()

# Load node list data
cout("Reading file...")
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
        body_nodes = sum(nl_id == n)
        body_mass = sum(m[nl_id == n])
        body_pos = pos[nl_id == n]
        rx = (body_pos[:,0].max() - body_pos[:,0].min())/2
        ry = (body_pos[:,1].max() - body_pos[:,1].min())/2
        rz = (body_pos[:,2].max() - body_pos[:,2].min())/2
        body_r = (rx + ry + rz)/3
        print "    List {:g}: {:.6g} kg in {} nodes ({:.4g} kg/node).".format(
            n, body_mass, body_nodes, body_mass/body_nodes)
        if args.print_radii:
            print "    List {:g}: {:g} m mean radius.".format(n, body_r)
except:
    raise StandardError("Could not read data from {}".format(
        args.filename))

