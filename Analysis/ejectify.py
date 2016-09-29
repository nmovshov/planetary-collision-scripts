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
        ejc, ind = ahelpers.ejectify_fnl(fnl, method=args.method)

        # Now we try to find original depth of ejecta particles
        orig_fnl = None
        origfile = os.path.basename(onefile)
        origfile = origfile[:origfile.find('-')] + '-00000-0.fnl'
        origfile = os.path.join(dirname, origfile)
        print "Looking for matching pre-impact file in {}.".format(dirname)
        try:
            orig_fnl = ahelpers.load_fnl(origfile)
            assert orig_fnl.nbNodes == fnl.nbNodes
        except (StandardError, AssertionError):
            try:
                orig_fnl = ahelpers.load_fnl(origfile + '.gz')
                assert orig_fnl.nbNodes == fnl.nbNodes
            except (StandardError, AssertionError):
                orig_fnl = None
                print ("Could not find valid pre-impact data;" +
                    " skipping calculation of initial depth.")

        R_ini = np.zeros(fnl.nbNodes)
        if orig_fnl is not None:
            for n in np.unique(orig_fnl.id):
                body = ahelpers.unpack_fnl(orig_fnl)[orig_fnl.id == n]
                pos = body[:,ahelpers.FNLMeta.x_col:ahelpers.FNLMeta.x_col+3]
                m = body[:,ahelpers.FNLMeta.m_col]
                X = [pos[:,0].dot(m), pos[:,1].dot(m), pos[:,2].dot(m)]/sum(m)
                pos = pos - X
                r = np.sqrt(pos[:,0]**2 + pos[:,1]**2 + pos[:,2]**2)
                R_ini[orig_fnl.id == n] = r/max(r)
                pass

        # Tack R_ini of ejecta to right of unpacked ejc
        ejc_arr = ahelpers.unpack_fnl(ejc)
        R_ini = R_ini[ind,None] # the None makes it n-by-1
        ejc_arr = np.hstack((ejc_arr,R_ini))

        # Delete useless hmax column
        ejc_arr = np.delete(ejc_arr, ahelpers.FNLMeta.hmax_col, 1)

        # Write header and save to file
        M_T = fnl.m.sum()
        M_LB = fnl.m.sum() - ejc.m.sum()
        head = ahelpers.fnl_header_ejecta.format(M_T,
                                                 fnl.nbNodes,
                                                 M_LB,
                                                 ejc.nbNodes)
        outname = os.path.join(dirname, 'ejecta_from_'+os.path.basename(onefile))
        format = 2*['%2d'] + (ahelpers.FNLMeta.nb_columns - 2)*['%12.5e']
        np.savetxt(outname, ejc_arr, header=head, fmt=format)
        #ahelpers.save_fnl(outname, ejc, head)
        print "Ejecta field saved to file {}".format(os.path.relpath(outname))
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

