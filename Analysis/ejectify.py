#!/soft/scipy_0.13.0/CentOS_6/bin/python
#---------------------------------------------------------------------------------
# ejectify - a utility for extracting the ejecta field from SPHERAL run output.
#
# Author: Naor Movshovitz (nmovshov at gee mail dot com)
#---------------------------------------------------------------------------------
import sys, os, glob
import numpy as np
import argparse
import ahelpers
from time import time
cout = sys.stdout.write

def _main():
    """Entry point when used as command line utility (recommended)."""

    # Parse command line arguments
    args = _PCL()

    # Ad hoc file-by-file treatment
    if os.path.isfile(args.filename):
        allfiles = [args.filename]
        dirname = os.path.dirname(os.path.abspath(args.filename))
    else:
        dirname = os.path.abspath(args.filename)
        allfiles = glob.glob(os.path.join(dirname, '*.fnl')) + \
                   glob.glob(os.path.join(dirname, '*.fnl.gz'))
        allfiles.sort()
    if len(allfiles) == 0:
        print "{} does not contain any valid fnl or fnl.gz files.".format(dirname)
        return
    
    ot = time()
    print
    for onefile in allfiles:
        # Load node list data
        cout("Reading file {}...".format(os.path.relpath(onefile)))
        try:
            fnl = ahelpers.load_fnl(onefile)
            cout("Done.\n")
            print "Found {2} kg in {1} node lists ({0} nodes total).".format(
                fnl.nbNodes, np.unique(fnl.id).size, sum(fnl.m))
            for n in np.unique(fnl.id):
                print "    List {:g}: {:.6g} kg in {} nodes ({:.4g} kg/node).".format(
                    n, sum(fnl.m[fnl.id == n]), sum(fnl.id == n),
                    sum(fnl.m[fnl.id == n])/sum(fnl.id == n))
        except StandardError:
            raise StandardError("Could not read data from {}".format(onefile))

        # Dispatch to the work method
        t = time()
        print "Ejectifying using algorithm {}...".format(args.method),
        sys.stdout.flush()
        tic = time()
        ejc = ahelpers.ejectify_fnl(fnl, method=args.method)
        print "Ejecta field saved to file {}".format('mitzi.ejecta')
        print "Elapsed time = {:g} sec.".format(time() - tic)
        print

    # Finish and exit
    print "All files done. Elapsed time = {:g} sec.".format(time() - ot)
    return

def _PCL():
    known_methods = ['jutzi', 'naor1']
    parser = argparse.ArgumentParser()
    parser.add_argument('filename',
        help="name of file containing node list data to ejectify")
    parser.add_argument('-m','--method',
        help="choice of algorithm for bound mass detection",
        choices=known_methods,
        default='naor1',
        action='store')
    #parser.add_argument('-q','--quiet',
    #    help="suppress progress output to stdout",
    #    action='store_true')
    parser.add_argument('-d','--delimiter',
        help="optional single-character delimiter for non FNL files",
        type=str,
        choices=[','],
        default=None)
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    _main()
    pass

