#------------------------------------------------------------------------------
#   Spheral Helpers - A module containing some convenience methods.
#
#------------------------------------------------------------------------------
import sys
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

    print 'Pickling', nl.label(), nl.name, '........'

    # Estimate memory usage and give user a chance to avoid crash
    nbFields = 11 # pos and vel count as 3 each
    bytesPerNode = 8*nbFields
    bytes = 2*bytesPerNode*nl.numNodes # times 2 b/c of temporary vars
    if bytes > 1e9:
        abort = raw_input('It looks like this will require a dangerous amount of memory.' +
                          '\nContinue anyway (%d bytes needed)? y/[n] '%bytes)
        if abort in ('n', 'no', 'N', 'No', 'NO', ''):
            return

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

    # Zip fields so that we have all fields for each node in the same tuple
    #  We do this so we can concatenate everything in a single reduction operation,
    #  to ensure that all fields in one record in the final list belong to the same
    #  node.
    localFields = zip(xloc, vloc, mloc, rloc, uloc, ploc, Tloc)

    # Do a SUM reduction on all ranks
    #  This works because the + operator for python lists is a concatenation!
    globalFields = mpi.allreduce(localFields, mpi.SUM)

    # Create a dictionary to hold field variables
    nlFieldDict = dict(x=[],   # position vector
                       v=[],   # velocity vector
                       m=[],   # mass
                       rho=[], # mass density
                       p=[],   # pressure
                       T=[],   # temperature
                       U=[],   # specific thermal energy
                      )

    # Loop over nodes to fill field values
    nbGlobalNodes = mpi.allreduce(nl.numInternalNodes, mpi.SUM)
    for k in range(nbGlobalNodes):
        nlFieldDict[  'x'].append((globalFields[k][0].x, globalFields[k][0].y, globalFields[k][0].z))
        nlFieldDict[  'v'].append((globalFields[k][1].x, globalFields[k][1].y, globalFields[k][1].z))
        nlFieldDict[  'm'].append( globalFields[k][2])
        nlFieldDict['rho'].append( globalFields[k][3])
        nlFieldDict[  'U'].append( globalFields[k][4])
        nlFieldDict[  'p'].append( globalFields[k][5])
        nlFieldDict[  'T'].append( globalFields[k][6])

    # Optionally, pickle the dict to a file
    if filename is not None:
        if isinstance(filename, str):
            with open(filename, 'wb') as fid:
                pickle.dump(nlFieldDict, fid)
        else:
            raise UserWarning('Dict NOT pickled to file because argument 2 is %s instead of %s' % 
                              (type(filename), type('x')))

    # And Bob's our uncle
    print 'Done.'
    return nlFieldDict
    # End function spickle_node_list

def pflatten_node_list(nl,filename):
    """Flatten physical field values from a node list to a rectangular ascii file.

    pflatten_node_list(nl,filename)

    See also: spickle_node_list
    """

    # And Bob's our uncle
    print 'alo world'
    return 0
    # End function pflatten_node_list

