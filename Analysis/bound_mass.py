#!/soft/scipy_0.13.0/CentOS_6/bin/python
#---------------------------------------------------------------------------------
# bound_mass - a utility for finding the largest gravitationally bound mass in SPH
# output data. 
#
# Author: Naor Movshovitz (nmovshov at gee mail dot com)
#---------------------------------------------------------------------------------
import sys, os, shutil
import numpy as np
import argparse
import ahelpers
from time import time
from numba import jit
cout = sys.stdout.write

def _main():
    """Entry point when used as command line utility (recommended)."""

    # Parse command line arguments
    args = _PCL()
    
    # Load node list data
    cout("Reading file...")
    try:
        fnl = ahelpers.load_fnl(args.filename)
        cout("Done.\n")
        pos = np.vstack((fnl.x, fnl.y, fnl.z)).T
        vel = np.vstack((fnl.vx, fnl.vy, fnl.vz)).T
        m   = fnl.m
        print "Found {2} kg in {1} node lists ({0} nodes total).".format(
            fnl.nbNodes, np.unique(fnl.id).size, sum(m))
        for n in np.unique(fnl.id):
            print "    List {:g}: {:.6g} kg in {} nodes ({:.4g} kg/node).".format(
                n, sum(m[fnl.id == n]), sum(fnl.id == n),
                sum(m[fnl.id == n])/sum(fnl.id == n))
    except StandardError:
        try:
            raw = np.loadtxt(args.filename, delimiter=args.delimiter)
            cout("Done.\n")
            pos = raw[:,0:3]
            vel = raw[:,3:6]
            m   = raw[:,6]
            print "Found {} particles totaling {} kg.".format(
                len(pos),sum(m))
        except:
            raise StandardError("Could not read data from {}".format(
                args.filename))

    # Dispatch to the work method
    print
    t = time()
    units = [1,1,1]
    for k in range(len(args.method)):
        print "Detecting bound mass using algorithm {}...".format(args.method[k])
        tic = time()
        [M_bound, ind_bound] = bound_mass(pos, vel, m,
                                          method=args.method[k],
                                          length_scale=args.length_scale,
                                          units=units)
        print "Found {:g} kg in {:g} particles; M_bound/M_tot = {:.4g}.".format(
            M_bound, sum(ind_bound), M_bound/sum(m))
        print "Elapsed time = {:g} sec.".format(time() - tic)
        print
    print "All done. Elapsed time = {:g} sec.".format(time() - t)

    # Exit
    return

def bound_mass(pos, vel, m, method='jutzi', length_scale=None, units=[1,1,1]):
    """Given cloud of particles return largest gravitationally bound mass.

    This function looks at a cloud of point masses with known positions and
    velocities and attempts to find the largest (by mass) subset that is
    gravitationally bound. There are different choices for the algorithm to use,
    and all algorithms are only approximations. The only sure way of finding the
    bound particles is to integrate the system in time for many gravitational time
    scales.

    The algorithm employed is chosen by the label supplied in the method argument.
    The choices implemented currently are:
      'kory1' - In the RF of the particle with lowest potential return the
                particles with negative total energy.
      'kory2' - In the RF with the most particles with negative total energy among
                all possible RFs centered on a particle return the particles with
                negative energy.
      'jutzi' - In the RF of the particle with lowest potential remove particles
                with positive total energy. Repeat until stable.
      'naor1' - Add particles bound to the particle with lowest potential. In the
                CM frame of this set, add particles bound to the set. Repeat until
                stable.
      'naor2' - In the RF centered on the CM of the largest spatially contiguous
                fragment return particles bound to that fragment.

    This function is a dispatcher - the work is carried out in sub functions.

    Parameters
    ----------
    pos : n-by-3 numeric array
        Particle positions.
    vel : n-by-3 numeric array
        Particle velocities.
    m : n-by-1 numeric array
        Particle masses.
    units : numeric positive 3-vector, optional
        Length, mass, and time units, in mks, of the particle coordinates.
    method : string
        Algorithm to use.
    length_scale : numeric, positive
        Override default length scale used in algorithm naor2

    Returns
    -------
    M_bound : real positive scalar
        Mass of the largest bound clump
    ind_bound : logical nparray
        Indices of bound particles
    """

    # Some minimal assertions (NOT bullet-proof filter!)
    pos = np.array(pos)
    vel = np.array(vel)
    m   = np.array(m)
    units = np.array(units)
    assert pos.ndim == 2 and pos.shape[1] == 3 and np.all(np.isreal(pos))
    assert vel.ndim == 2 and vel.shape[1] == 3 and np.all(np.isreal(vel))
    assert m.ndim == 1 and np.all(np.isreal(m)) and np.all(m > 0)
    assert units.ndim == 1 and len(units) == 3 and np.all(units > 0)
    assert len(pos) == len(vel) == len(m)
    assert method in ['kory1', 'kory2', 'jutzi', 'naor1', 'naor2']
    assert np.size(length_scale) == 1 and np.isreal(length_scale)

    # Deal with units
    bigG = 6.67384e-11*units[0]**(-3)*units[1]*units[2]**2
    length_scale = length_scale*units[0]

    # Dispatch to sub functions by method
    if   method == 'kory1':
        (M_bound, ind_bound) = _bm_kory1(pos, vel, m, bigG)
        pass
    elif method == 'kory2':
        (M_bound, ind_bound) = _bm_kory2(pos, vel, m, bigG)
        pass
    elif method == 'jutzi':
        (M_bound, ind_bound) = _bm_jutzi(pos, vel, m, bigG)
        pass
    elif method == 'naor1':
        (M_bound, ind_bound) = _bm_naor1(pos, vel, m, bigG)
        pass
    elif method == 'naor2':
        (M_bound, ind_bound) = _bm_naor2(pos, vel, m, bigG, length_scale)
        pass
    else:
        sys.exit("Unknown method") # this can't really happen
        pass

    return (M_bound, ind_bound)

def _bm_kory1(pos, vel, m, bigG):
    """In RF of lowest potential node return nodes with negative energy."""
    U = bigG*_potential(pos[:,0], pos[:,1], pos[:,2], m);
    ind = np.argmin(U)
    VCM = vel[ind,:]
    ind_bound = np.array(len(m)*[False])
    for j in range(len(m)):
        V = vel[j,:] - VCM
        K = 0.5*(V[0]*V[0] + V[1]*V[1] + V[2]*V[2])
        if K + U[j] < 0:
            ind_bound[j] = True
            pass
        pass
    return (sum(m[ind_bound]), ind_bound)

def _bm_kory2(pos, vel, m, bigG):
    """Use RF with most bound nodes among all possible RFs centered on a node."""
    U = bigG*_potential(pos[:,0], pos[:,1], pos[:,2], m);
    max_M = -np.inf
    ind_bound = np.array(len(m)*[False])
    K = np.zeros(len(m))
    for j in range(len(m)):
        VCM = vel[j]
        for k in range(len(m)):
            V = vel[j,:] - VCM
            K[k] = 0.5*(V[0]*V[0] + V[1]*V[1] + V[2]*V[2])
            pass
        mask = K + U < 0
        if sum(m[mask]) > max_M:
            ind_bound = mask
            pass
        pass
    return (sum(m[ind_bound]), ind_bound)

def _bm_jutzi(pos, vel, m, bigG):
    """In RF of lowest potential remove nodes with positive energy and repeat."""
    bU = bigG*_potential(pos[:,0], pos[:,1], pos[:,2], m)
    ind = np.argmin(bU)
    VCM = vel[ind,:]
    ind_bound = np.array(len(m)*[True])
    nbb = -1
    while nbb != sum(ind_bound):
        nbb = sum(ind_bound)
        bU = bigG*_potential(pos[:,0], pos[:,1], pos[:,2], m, ind_bound)
        for j in range(len(m)):
            V = vel[j,:] - VCM
            K = 0.5*(V[0]*V[0] + V[1]*V[1] + V[2]*V[2])
            if K + bU[j] > 0.0:
                ind_bound[j] = False
                pass
            pass
        pass
    return (sum(m[ind_bound]), ind_bound)

def _bm_naor1(pos, vel, m, bigG):
    """Add nodes bound to CM of bound nodes until stable. Seed with lowest U."""
    ind_bound = np.array(len(m)*[False])
    U = bigG*_potential(pos[:,0], pos[:,1], pos[:,2], m)
    ind = np.argmin(U)
    ind_bound[ind] = True
    nbb = -1
    m3 = np.tile(m,(3,1)).T
    while nbb != sum(ind_bound):
        nbb = sum(ind_bound)
        M = sum(m[ind_bound])
        cmpos = np.sum(m3[ind_bound,:]*pos[ind_bound,:], 0)/M
        cmvel = np.sum(m3[ind_bound,:]*vel[ind_bound,:], 0)/M
        for j in range(len(m)):
            dR = pos[j,:] - cmpos
            dr = np.sqrt(np.dot(dR, dR)) + np.spacing(1)
            U = -bigG*M/dr
            V = vel[j,:] - cmvel
            K = 0.5*(V[0]*V[0] + V[1]*V[1] + V[2]*V[2])
            if K + U < 0.0:
                ind_bound[j] = True
            else:
                ind_bound[j] = False
            pass
        pass
    pass
    return (sum(m[ind_bound]), ind_bound)

def _bm_naor2(pos, vel, m, bigG, length_scale):
    """Add nodes bound to CM of largest spatially contiguous clump."""
    
    # Deal with length scale
    if length_scale == 0.0:
        length_scale = max((pos.max(0) - pos.min(0))/m.size)

    # First, find the largest clump based on euclidean proximity. Note that this
    # is the time consuming part of the process, using an n^2 algorithm for
    # neighbor finding.
    labels = fast_clumps(pos.tolist(), length_scale)
    labels = np.array(labels)
    c_labels = np.unique(labels)
    c_masses = [sum(m[labels == c_labels[k]]) for k in range(len(c_labels))]
    c_major_label = c_labels[np.argmax(c_masses)]
    
    # Ok. All nodes labeled c_major_label belong to the biggest geometrical
    # clump. Let's find the center of mass of this clump.
    cmask = labels == c_major_label
    POS = pos[cmask]
    VEL = vel[cmask]
    M = m[cmask]
    R_com = np.dot(POS.T, M)/sum(M)
    V_com = np.dot(VEL.T, M)/sum(M)

    # Ok. Now pretend each node outside the largest clump feels a point-mass
    # potential from the clump, and test its energy in the clump frame.
    M_bound = M.sum()
    ind_bound = cmask
    GM = bigG*M_bound
    for k in range(len(pos)):
        if ind_bound[k]:
            continue
        U = -GM/np.sqrt(np.dot(pos[k], pos[k]))
        K = np.dot(vel[k], vel[k])
        if U + K < 0:
            M_bound += m[k]
            ind_bound[k] = True
            pass
        pass

    # That's it.
    return (M_bound, ind_bound)

@jit
def fast_clumps(pos, L):
    """Partition a cloud of point masses into distinct clumps based on proximity.

    This is the optimized version of nr3.eclazz functions, which partitions any
    set (tree) into equivalence classes (connected components). See nr3.eclazz for
    details of the algorithm. This function is a specialized version for the case
    where the set is a cloud of particles with (x,y,z) coordinates and the
    connectivity test is a simple euclidean distance threshold.

    Parameters
    ----------
    pos : numpy.ndarray or list
      An n-by-3 array of coordinates. A list of lists seems to work best.
    L : float, positive
      The distance threshold for the proximity test.

    Returns
    -------
    labels : list of ints
        A vector of integer labels. labels[k] is the clump that pos[k] belongs to.
        There are len(np.unique(labels)) such clumps.
    """
    
    ## Some minimal assertions
    #assert isinstance(pos, (list, np.ndarray))
    #assert np.isreal(L) and np.isscalar(L) and L >= 0
    #if type(pos) is np.ndarray:
    #    assert(pos.ndim == 2 and pos.shape[1] == 3)
    #    pass

    # Prepare
    labels = [-1]*len(pos)
    L2 = L**2
    
    # This is the nr3 algorithm, specialized to a euclidean distance test
    labels[0] = 0;
    for j in range(1, len(pos)):
        labels[j] = j
        pj = pos[j]
        for k in range(0,j):
            labels[k] = labels[labels[k]]
            pk = pos[k]
            if (pj[0] - pk[0])**2 + (pj[1] - pk[1])**2 + (pj[2] - pk[2])**2 < L2:
                labels[labels[labels[k]]] = j
                pass
            pass
        pass
    for j in range(len(pos)):
        labels[j] = labels[labels[j]]
        pass
    
    # That's it
    return labels

def _PCL():
    known_methods = ['kory1', 'kory2', 'jutzi', 'naor1', 'naor2']
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help="name of file containing node list data")
    parser.add_argument('-m','--method',
        help="choice of algorithm; may be specified multiple times",
        choices=known_methods + ['all'],
        default=[],
        action='append')
    #parser.add_argument('-q','--quiet',
    #    help="suppress progress output to stdout",
    #    action='store_true')
    parser.add_argument('-I','--max-iter',
        help="max number of iterations in iterative methods",
        type=int,
        default=10)
    parser.add_argument('-d','--delimiter',
        help="optional single-character delimiter for non FNL files",
        type=str,
        choices=[','],
        default=None)
    parser.add_argument('-L','--length-scale',
        help="length scale in meters for proximity test",
        type=float,
        default=0.0)
    args = parser.parse_args()
    if 'all' in args.method:
        args.method = known_methods
    elif args.method == []:
        args.method = ['naor1', 'jutzi']
    else:
        pass
    return args

def _potential(x, y, z, m, mask=None):
    if mask is None:
        mask = np.array(len(x)*[True])
    return _c_potential(x, y, z, m, mask)

@jit
def _c_potential(x, y, z, m, mask):
    U = np.zeros(x.shape)
    for j in range(len(x)):
        if mask[j]:
            for k in range(j):
                if mask[k]:
                    dx = x[j] - x[k]
                    dy = y[j] - y[k]
                    dz = z[j] - z[k]
                    dr = np.sqrt(dx*dx + dy*dy + dz*dz)
                    U[j] = U[j] - m[k]/dr
                    U[k] = U[k] - m[j]/dr
                    pass
                pass
            pass
        pass
    return U

if __name__ == "__main__":
    _main()
    pass

