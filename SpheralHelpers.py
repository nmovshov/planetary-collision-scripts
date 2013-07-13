#------------------------------------------------------------------------------
#   Spheral Helpers - A module containing some convenience methods
#
#       * pickle_node_list - saves the field variables of a node list
#------------------------------------------------------------------------------
import sys, mpi
import cPickle as pickle
from SolidSpheral3d import *

def pickle_node_list(nl,filename=None):
    """Flatten a node list, saving the physical field variables in a dict.

    (This is not a true pickler class.)

    pickle_node_list(nl,filename) extracts field variables from nl which must
    be a valid node list, packs them in a dict, and attempts to pickle it to a
    file filename. The file will be overwritten if it exists. Dict is returned
    to caller.
    """

    print 'Pickling', nl.label(), nl.name, '...'

    # get references to field variables stored in node list
    xref = nl.positions()
    vref = nl.velocity()
    mref = nl.mass()
    rref = nl.massDensity()
    uref = nl.specificThermalEnergy()

    # pressure and temperature are stored in the eos object
    eos = nl.equationOfState()
    pref = ScalarField('pref',nl)
    Tref = ScalarField('Tref',nl)
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


