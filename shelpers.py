#------------------------------------------------------------------------------
#   Spheral Helpers - A collection of some convenience functions for reuse in
#                     the planetary collision scripts.
#
# Author: nmovshov@gmail.com
#------------------------------------------------------------------------------
import sys, os
import warnings
import mpi # Mike's simplified mpi wrapper
import cPickle as pickle
import SolidSpheral3d as sph

def spickle_node_list(nl,filename=None):
    """Pack physical field variables from a node list in a dict and pickle.

    (Note: This is not a true pickler class.)

    spickle_node_list(nl,filename) extracts field variables from all nodes of nl,
    which must be a valid node list, and packs them in a dict that is returned
    to the caller. If the optional argument filename is a string then dict will
    also be pickled to a file of that name. The file will be overwritten if it
    exists.

    The s in spickle is for 'serial', a reminder that this method collects ALL
    nodes of the node list (not just local nodes) in a single process. Thus this
    method is mainly useful for interactive work with small node lists. It is the
    user's responsibility to make sure her process has enough memory to hold the
    returned dict.

    See also: pflatten_node_list
    """

    # Make sure we are not wasting our time.
    assert isinstance(nl,(sph.Spheral.NodeSpace.FluidNodeList3d,
                          sph.Spheral.SolidMaterial.SolidNodeList3d)
                     ), "argument 1 must be a node list"
    
    # Estimate memory usage and give user a chance to avoid a crash.
    # TODO fix this block
    if mpi.rank == 0:
        nbFields = 11 # pos and vel count as 3 each
        bytesPerNode = 8*nbFields
        bytes = 2*bytesPerNode*nl.numNodes # times 2 b/c of temporary vars
        if bytes > 1e9:
            msg = "It looks like this will require a dangerous amount of memory." + \
                  "\nContinue anyway (%d bytes needed)? y/[n] " % bytes
            abort = raw_input(msg)
            if abort in ('n', 'no', 'N', 'No', 'NO', ''):
                return
            pass
        pass
    mpi.barrier()

    # Start collecting data.
    print 'Pickling', nl.label(), nl.name, '........'

    # Get values of field variables stored in internal nodes.
    xloc = nl.positions().internalValues()
    vloc = nl.velocity().internalValues()
    mloc = nl.mass().internalValues()
    rloc = nl.massDensity().internalValues()
    uloc = nl.specificThermalEnergy().internalValues()
    #(pressure and temperature are stored in the eos object and accessed differently)
    eos = nl.equationOfState()
    ploc = sph.ScalarField('ploc',nl)
    Tloc = sph.ScalarField('loc',nl)
    rref = nl.massDensity()
    uref = nl.specificThermalEnergy()
    eos.setPressure(ploc,rref,uref)
    eos.setTemperature(Tloc,rref,uref)

    # Zip fields so that we have all fields for each node in the same tuple.
    #  We do this so we can concatenate everything in a single reduction operation,
    #  to ensure that all fields in one record in the final list belong to the same
    #  node.
    localFields = zip(xloc, vloc, mloc, rloc, uloc, ploc, Tloc)

    # Do a SUM reduction on all ranks.
    #  This works because the + operator for python lists is a concatenation!
    globalFields = mpi.allreduce(localFields, mpi.SUM)

    # Create a dictionary to hold field variables.
    nlFieldDict = dict(x=[],   # position vector
                       v=[],   # velocity vector
                       m=[],   # mass
                       rho=[], # mass density
                       p=[],   # pressure
                       T=[],   # temperature
                       U=[],   # specific thermal energy
                      )

    # Loop over nodes to fill field values.
    nbGlobalNodes = mpi.allreduce(nl.numInternalNodes, mpi.SUM)
    for k in range(nbGlobalNodes):
        nlFieldDict[  'x'].append((globalFields[k][0].x, globalFields[k][0].y, globalFields[k][0].z))
        nlFieldDict[  'v'].append((globalFields[k][1].x, globalFields[k][1].y, globalFields[k][1].z))
        nlFieldDict[  'm'].append( globalFields[k][2])
        nlFieldDict['rho'].append( globalFields[k][3])
        nlFieldDict[  'U'].append( globalFields[k][4])
        nlFieldDict[  'p'].append( globalFields[k][5])
        nlFieldDict[  'T'].append( globalFields[k][6])

    # Optionally, pickle the dict to a file.
    if mpi.rank == 0:
        if filename is not None:
            if isinstance(filename, str):
                with open(filename, 'wb') as fid:
                    pickle.dump(nlFieldDict, fid)
                    pass
                pass
            else:
                msg = "Dict NOT pickled to file because " + \
                      "argument 2 is {} instead of {}".format(type(filename), type('x'))
                sys.stderr.write(msg+'\n')
                pass
            pass
        pass
        
    # And Bob's our uncle.
    print "Done."
    return nlFieldDict
    # End function spickle_node_list


def pflatten_node_list(nl,filename,do_header=True,nl_id=0):
    """Flatten physical field values from a node list to a rectangular ascii file.

    pflatten_node_list(nl,filename) extracts field variables from all nodes of nl,
    which must be a valid node list, and writes them as a rectangular table into
    the text file filename. (A short header is also written, using the # comment
    character so the resulting file can be easily read with, e.g., numpy.loadtext.)
    The file will be overwritten if it exists.

    pflatten_node_list(...,do_header=False) omits the header and appends the flattened
    nl to the end of the files.

    pflatten_node_list(...,nl_id=id) inserts the integer id in the first column
    of every node (row) in the node list. This can be used when appending multiple
    lists to the same file, providing a convenient way to distinguish nodes from
    different lists when the file is later read.

    The format of the output table is (one line per node):
      x y z vx vy vz m rho p T U

    The p in pflatten is for 'parallel', a reminder that the all nodes will be
    processed in their local rank, without ever being communicated or collected
    in a single process. Each mpi rank will wait its turn to access the output file,
    so the writing is in fact serial, but avoids bandwidth and memory waste and
    is thus suitable for large node lists from high-res runs.

    See also: spickle_node_list
    """

    # Make sure we are not wasting our time.
    assert isinstance(nl,(sph.Spheral.NodeSpace.FluidNodeList3d,
                          sph.Spheral.SolidMaterial.SolidNodeList3d)
                     ), "argument 1 must be a node list"
    assert isinstance(filename, str), "argument 2 must be a simple string"
    assert isinstance(do_header, bool), "true or false"
    assert isinstance(nl_id, int), "int only idents"
    assert not isinstance(nl_id, bool), "int only idents"

    # Write the header
    if do_header:
        if mpi.rank == 0:
            fid = open(filename,'w')
            header = "#Place holder for real header\n"
            fid.write(header)
            fid.close()
            pass
        pass
     
    # Start collecting data.
    print 'Flattening', nl.label(), nl.name, '........'
    
    # Get values of field variables stored in internal nodes
    xloc = nl.positions().internalValues()
    vloc = nl.velocity().internalValues()
    mloc = nl.mass().internalValues()
    rloc = nl.massDensity().internalValues()
    uloc = nl.specificThermalEnergy().internalValues()
    #(pressure and temperature are stored in the eos object and accessed differently)
    eos = nl.equationOfState()
    ploc = sph.ScalarField('ploc',nl)
    Tloc = sph.ScalarField('loc',nl)
    rref = nl.massDensity()
    uref = nl.specificThermalEnergy()
    eos.setPressure(ploc,rref,uref)
    eos.setTemperature(Tloc,rref,uref)

    # Procs take turns writing internal node values to file.
    for proc in range(mpi.procs):
        if proc == mpi.rank:
            fid = open(filename,'a')
            for nk in range(nl.numInternalNodes):
                line  = "{:2d}  ".format(nl_id)
                line += "{0.x:+12.5e}  {0.y:+12.5e}  {0.z:+12.5e}  ".format(xloc[nk])
                line += "{0.x:+12.5e}  {0.y:+12.5e}  {0.z:+12.5e}  ".format(vloc[nk])
                line += "{0:+12.5e}  ".format(mloc[nk])
                line += "{0:+12.5e}  ".format(rloc[nk])
                line += "{0:+12.5e}  ".format(uloc[nk])
                line += "{0:+12.5e}  ".format(ploc[nk])
                line += "{0:+12.5e}  ".format(Tloc[nk])
                line += "\n"
                fid.write(line)
                pass
            fid.close()
            pass
        mpi.barrier()
        pass
     
    # And Bob's our uncle
    print "Done."
    # End function pflatten_node_list
     

