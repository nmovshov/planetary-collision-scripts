#------------------------------------------------------------------------------
#   Spheral Helpers - A module containing some convenience methods.
#
#------------------------------------------------------------------------------
import sys, mpi
import cPickle as pickle
import SolidSpheral3d as sph

def spickle_node_list(nl,filename=None):
    """Flatten physical field variables from a node list and return in a dict.

    (Note: This is not a real pickler class.)

    spickle_node_list(nl,filename) extracts field variables from ALL nodes of nl,
    which must be a valid node list, and packs them in a dict that is returned
    to the caller. If the optional argument filename is a string the dict will
    also be pickled to a file of that name in the current directory. The file
    will be overwritten if it exists.

    The s in spickle is for 'serial' because this method collects ALL nodes of
    the node list (not just local nodes) in a single process. Thus this method
    is useful mostly for interactive work with small node lists. It is the user's
    responsibility to make sure her process has enough memory to hold the returned
    dict.

    See also: ppickle_node_list
    """

    print 'Pickling', nl.label(), nl.name, '...'

    # estimate memory usage and give user a chance to avoid crash
    nbFields = 11 # pos and vel count as 3 each
    bytesPerNode = 8*nbFields
    bytes = bytesPerNode*nl.numNodes
    if bytes > 2e2:
        abort = raw_input('It looks like this will require a dangerous amount of memory.' +
                          '\nContinue anyway (%d bytes needed)? y/[n] '%bytes)
        if abort in ('n', 'no', 'N', 'No', 'NO', ''):
            return

    # get references to field variables stored in node list
    xref = nl.positions()
    vref = nl.velocity()
    mref = nl.mass()
    rref = nl.massDensity()
    uref = nl.specificThermalEnergy()

    # pressure and temperature are stored in the eos object
    eos = nl.equationOfState()
    pref = sph.ScalarField('pref',nl)
    Tref = sph.ScalarField('Tref',nl)
    eos.setPressure(pref,rref,uref)
    eos.setTemperature(Tref,rref,uref)

    # create a dictionary to hold field variables
    nlfields = dict(x=[],   # position vector
                    v=[],   # velocity vector
                    m=[],   # mass
                    rho=[], # mass density
                    P=[],   # pressure
                    T=[],   # temperature
                    U=[],   # specific thermal energy
                    )

    # loop over nodes to fill field values
    for k in range(nl.numInternalNodes):
        nlfields['x'].append((xref[k].x, xref[k].y, xref[k].z))
        nlfields['v'].append((vref[k].x, vref[k].y, vref[k].z))
        nlfields['m'].append(mref[k])
        nlfields['rho'].append(rref[k])
        nlfields['P'].append(pref[k])
        nlfields['T'].append(Tref[k])
        nlfields['U'].append(uref[k])

    print 'Done.'
    return nlfields


