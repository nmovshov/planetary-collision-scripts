#-------------------------------------------------------------------------------
#   Spheral Helpers - A collection of some convenience functions for reuse in
#                     the planetary collision scripts.
#
# Author: nmovshov at gmail dot com
#-------------------------------------------------------------------------------
import sys, os
import mpi # Mike's simplified mpi wrapper
import cPickle as pickle
import numpy as np
import scipy as sp
import SolidSpheral3d as sph

def hydrostaticize_two_layer_planet(inner, outer, G=6.674e-11):
    """Modify densities in node generators to approximate hydrostatic equilibrium.

    Assuming a barely compressible, two-layer planet, a pressure profile in
    hydrostatic equilibrium can be found by integrating the hydrostatic equation
    with constant density. The equation of state can then be inverted to provide
    a density profile consistent with this pressure profile. Although the 
    resulting pressure/density state is not strictly self consistent, it may be
    used as a good approximation for small planets that are not expected to be
    highly compressed.

    This function takes in two node generators of the hcp class, and modifies the
    density and mass of nodes in each to match a hydrostatic state. To invert the
    equation of state this function uses the scipy solver.
    """

    # Make sure we are not wasting our time.
    import PlanetNodeGenerators as PNG
    assert isinstance(inner, PNG.HexagonalClosePacking), "must be HCP generator"
    assert isinstance(outer, PNG.HexagonalClosePacking), "must be HCP generator"

    # Setup local variables
    R = outer.rMax
    rc = inner.rMax
    rhoc = inner.rho[0]
    rhom = outer.rho[0]
    assert 0 < rc < R
    assert rhom < rhoc
    r_inner = np.hypot(inner.x, np.hypot(inner.y, inner.z))
    r_outer = np.hypot(outer.x, np.hypot(outer.y, outer.z))

    # Step one - calculate pressure profile
    c2 = 4*np.pi/3*G*(0.5*rhom**2*R**2 - rhom*(rhoc - rhom)*rc**3/R)
    c1 = 4*np.pi/3*G*(0.5*rhoc**2 - 1.5*rhom**2 + rhoc*rhom)*rc**2 + c2
    p_inner = np.ones(r_inner.size)*np.NaN
    p_outer = np.ones(r_outer.size)*np.NaN
    for k in range(r_inner.size):
        p_inner[k] = c1 - 4*np.pi/3*G*0.5*rhoc**2*r_inner[k]**2
    for k in range(r_outer.size):
        p_outer[k] = c2 - 4*np.pi/3*G*(0.5*rhom**2*r_outer[k]**2 - 
                                       rhom*(rhoc - rhom)*rc**3/r_outer[k])
    assert np.all(np.isfinite(p_inner))
    assert np.all(np.isfinite(p_outer))

    # Step two - use root finding to invert eos and get a density
    from scipy.optimize import brentq
    def f(x,p,eos):
        return p - eos.pressure(x,0)
    for k in range(r_inner.size):
        eos = inner.EOS
        inner.rho[k] = brentq(f, eos.referenceDensity/2, eos.referenceDensity*2, 
                              args=(p_inner[k],eos))
        inner.m[k] = inner.rho[k]*inner.V[k]
    for k in range(r_outer.size):
        eos = outer.EOS
        outer.rho[k] = brentq(f, eos.referenceDensity/2, eos.referenceDensity*2, 
                              args=(p_outer[k],eos))
        outer.m[k] = outer.rho[k]*outer.V[k]
    assert np.all(np.isfinite(inner.rho))
    assert np.all(np.isfinite(outer.rho))

    # And Bob's our uncle
    return 
    # End function hydrostaticize_two_layer_planet
    

def construct_eos_for_material(material_tag,units,etamin=0.94,etamax=100.0):
    """Return a spheral EOS object for a material identified by tag.

    construct_eos_for_material(mtag,units) calls the appropriate spheral eos
    constructor for the material identified by mtag, which must be one of the keys
    defined in the global shelpers.material_dictionary. This dictionary also 
    includes additional arguments to be passed to the constructor, when necessary.

    The etamin and etamax optional arguments have slightly different meaning 
    depending on which EOS constructor is actually used. Currently implemented
    constructors are:
      Tillotson : the value of etamin is passed to the etamin_solid parameter of
                  the constructor. This is used to limit tensional pressure when
                  the material is no longer solid. (Note that the spheral 
                  constructor also has an etamin parameter, which is used to 
                  prevent underflows in the pressure computation.)
      ANEOS : Not yet implemented.

    All uss runs should use this method to create equations of state, instead of
    calling the spheral constructors directly, in order to allow automatic record
    keeping of what material was used in a given run. This also allows reusing 
    "pre cooked" node lists in new runs.

    The file <uss>/MATERIALS.md should contain a table of available material tags.

    See also: material_dictionary
    """

    # Make sure we are not wasting our time.
    assert isinstance(material_tag,str)
    assert material_tag.lower() in material_dictionary.keys()
    assert isinstance(units,sph.PhysicalConstants)
    assert isinstance(etamin,float)
    assert isinstance(etamax,float)

    # Build eos using our internal dictionary
    mat_dict = material_dictionary[material_tag.lower()]
    eos_constructor = mat_dict['eos_constructor']
    eos_arguments = mat_dict['eos_arguments']
    eos = None

    if mat_dict['eos_type'] == 'Tillotson':
        eos = eos_constructor(eos_arguments['materialName'],
                              1e-20, 1e20, units,
                              etamin_solid=etamin)
        eos.uid = mat_dict['eos_id']
        pass
    else:
        print "EOS type {} not yet implemented".format(mat_dict['eos_type'])
        pass

    # And Bob's our uncle
    return eos
    # End function construct_eos_for_material


def spickle_node_list(nl,filename=None,silent=False):
    """Pack physical field variables from a node list in a dict and pickle.

    (Note: This is not a true pickler class.)

    spickle_node_list(nl,filename) extracts field variables from all nodes of nl,
    which must be a valid node list, and packs them in a dict that is returned
    to the caller. If the optional argument filename is a string then dict will
    also be pickled to a file of that name. The file will be overwritten if it
    exists.

    The s in spickle is for 'serial', a reminder that this method collects all
    nodes of the node list (from all ranks) in a single process. Thus this method
    is mainly useful for interactive work with small node lists. It is the user's
    responsibility to make sure her process has enough memory to hold the returned
    dict.

    See also: pflatten_node_list
    """

    # Make sure we are not wasting our time.
    assert isinstance(nl,(sph.Spheral.NodeSpace.FluidNodeList3d,
                          sph.Spheral.SolidMaterial.SolidNodeList3d)
                      ), "argument 1 must be a node list"
    assert isinstance(silent, bool), "true or false"
    
    # Start collecting data.
    if not silent:
        sys.stdout.write('Pickling ' +  nl.label() + ' ' + nl.name + '........')

    # Get values of field variables stored in internal nodes.
    xloc = nl.positions().internalValues()
    vloc = nl.velocity().internalValues()
    mloc = nl.mass().internalValues()
    rloc = nl.massDensity().internalValues()
    uloc = nl.specificThermalEnergy().internalValues()
    Hloc = nl.Hfield().internalValues()
    #(pressure and temperature are stored in the eos object.)
    eos = nl.equationOfState()
    ploc = sph.ScalarField('ploc',nl)
    Tloc = sph.ScalarField('loc',nl)
    rref = nl.massDensity()
    uref = nl.specificThermalEnergy()
    eos.setPressure(ploc,rref,uref)
    eos.setTemperature(Tloc,rref,uref)

    # Zip fields so that we have all fields for each node in the same tuple.
    #  We do this so we can concatenate everything in a single reduction operation,
    #  to ensure that all fields in one record in the final list belong to the 
    #  same node.
    localFields = zip(xloc, vloc, mloc, rloc, uloc, ploc, Tloc, Hloc)

    # Do a SUM reduction on all ranks.
    #  This works because the + operator for python lists is a concatenation!
    globalFields = mpi.allreduce(localFields, mpi.SUM)

    # Create a dictionary to hold field variables.
    nlFieldDict = dict(name=nl.name,
                       x=[],   # position vector
                       v=[],   # velocity vector
                       m=[],   # mass
                       rho=[], # mass density
                       p=[],   # pressure
                       h=[],   # smoothing ellipsoid axes
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
        nlFieldDict[  'h'].append((globalFields[k][7].Inverse().eigenValues().x,
                                   globalFields[k][7].Inverse().eigenValues().y,
                                   globalFields[k][7].Inverse().eigenValues().z,
                                   ))

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
    if not silent:
        print "Done."
    return nlFieldDict
    # End function spickle_node_list


def pflatten_node_list(nl,filename,do_header=True,nl_id=0,silent=False):
    """Flatten physical field values from a node list to a rectangular ascii file.

    pflatten_node_list(nl,filename) extracts field variables from all nodes of nl,
    which must be a valid node list, and writes them as a rectangular table into
    the text file filename. (A short header is also written, using the # comment
    character so the resulting file can be easily read with numpy.loadtext.) The
    file will be overwritten if it exists. If filename has the .gz extension it
    will be compressed using gzip.

    pflatten_node_list(...,do_header=False) omits the header and appends the 
    flattened nl to the end of the file if one exists.

    pflatten_node_list(...,nl_id=id) places the integer id in the first column
    of every node (row) in the node list. This can be used when appending multiple
    lists to the same file, providing a convenient way to distinguish nodes from
    different lists when the file is later read. The default id (for single node
    list files) is 0.

    The format of the output table is (one line per node):
      id eos_id x y z vx vy vz m rho p T U hmin hmax

    The p in pflatten is for 'parallel', a reminder that all nodes will be
    processed in their local rank, without ever being communicated or collected
    in a single process. Each mpi rank will wait its turn to access the output 
    file, so the writing is in fact serial, but avoids bandwidth and memory waste
    and is thus suitable for large node lists from high-res runs.

    See also: spickle_node_list
    """

    # Make sure we are not wasting our time.
    assert isinstance(nl,(sph.Spheral.NodeSpace.FluidNodeList3d,
                          sph.Spheral.SolidMaterial.SolidNodeList3d)
                      ), "argument 1 must be a node list"
    assert isinstance(filename, str), "argument 2 must be a simple string"
    assert isinstance(do_header, bool), "true or false"
    assert isinstance(silent, bool), "true or false"
    assert isinstance(nl_id, int), "int only idents"
    assert not isinstance(nl_id, bool), "int only idents"

    # Determine if file should be compressed.
    if os.path.splitext(filename)[1] == '.gz':
        import gzip
        open = gzip.open
    else:
        import __builtin__
        open = __builtin__.open

    # Write the header.
    if do_header:
        nbGlobalNodes = mpi.allreduce(nl.numInternalNodes, mpi.SUM)
        header = header_template.format(nbGlobalNodes)
        if mpi.rank == 0:
            fid = open(filename,'w')
            fid.write(header)
            fid.close()
            pass
        pass
     
    # Start collecting data.
    if not silent:
        sys.stdout.write('Flattening ' + nl.label() + ' ' + nl.name + '........')
    
    # Get values of field variables stored in internal nodes.
    xloc = nl.positions().internalValues()
    vloc = nl.velocity().internalValues()
    mloc = nl.mass().internalValues()
    rloc = nl.massDensity().internalValues()
    uloc = nl.specificThermalEnergy().internalValues()
    Hloc = nl.Hfield().internalValues()
    #(pressure and temperature are stored in the eos object.)
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
                line += "{:2d}  ".format(getattr(nl,'eos_id',-1))
                line += "{0.x:+12.5e}  {0.y:+12.5e}  {0.z:+12.5e}  ".format(xloc[nk])
                line += "{0.x:+12.5e}  {0.y:+12.5e}  {0.z:+12.5e}  ".format(vloc[nk])
                line += "{0:+12.5e}  ".format(mloc[nk])
                line += "{0:+12.5e}  ".format(rloc[nk])
                line += "{0:+12.5e}  ".format(ploc[nk])
                line += "{0:+12.5e}  ".format(Tloc[nk])
                line += "{0:+12.5e}  ".format(uloc[nk])
                line += "{0:+12.5e}  ".format(Hloc[nk].Inverse().eigenValues().minElement())
                line += "{0:+12.5e}  ".format(Hloc[nk].Inverse().eigenValues().maxElement())
                line += "\n"
                fid.write(line)
                pass
            fid.close()
            pass
        mpi.barrier()
        pass
     
    # And Bob's our uncle.
    if not silent:
        print "Done."
    # End function pflatten_node_list


def pflatten_node_list_list(nls,filename,do_header=True,silent=False):
    """Flatten a list of node lists to a rectangular ascii file.

    pflatten_node_list_list(nls,filename) writes meta data about the node lists
    in nls, which must be either a list or a tuple of valid node lists, to a 
    header of the file filename, and then calls pflatten_node_list(nl,filename)
    for each nl in nls.

    pflatten_node_list_list(...,do_header=False) omits the header.

    See also: pflatten_node_list
    """

    # Make sure we are not wasting our time.
    assert isinstance(nls,(list,tuple)), "argument 1 must be a list or tuple"
    assert isinstance(filename, str), "argument 2 must be a simple string"
    assert isinstance(do_header, bool), "true or false"
    assert isinstance(silent, bool), "true or false"
    for nl in nls:
        assert isinstance(nl,(sph.Spheral.NodeSpace.FluidNodeList3d,
                                  sph.Spheral.SolidMaterial.SolidNodeList3d)
                         ), "argument 1 must contain node lists"

    # Determine if file should be compressed.
    if os.path.splitext(filename)[1] == '.gz':
        import gzip
        open = gzip.open
    else:
        import __builtin__
        open = __builtin__.open

    # Write the header.
    if do_header:
        nbGlobalNodes = 0
        for nl in nls:
            nbGlobalNodes += mpi.allreduce(nl.numInternalNodes, mpi.SUM)
        header = header_template.format(nbGlobalNodes)
        if mpi.rank == 0:
            fid = open(filename,'w')
            fid.write(header)
            fid.close()
            pass
        pass

    # Send contents of nls to be flattened.
    for k in range(len(nls)):
        pflatten_node_list(nls[k],filename,do_header=False,nl_id=k,silent=silent)
        pass

    # And Bob's our uncle.
    if not silent:
        print "Done."
    # End function pflatten_node_list_list


global nb_fnl_columns
nb_fnl_columns = 15

global header_template
header_template = """\
################################################################################
# This file contains output data from a Spheral++ simulation, including all 
# field variables as well as some diagnostic data and node meta data. This
# file should contain {} data lines, one per SPH node used in the simulation.
# Line order is not significant and is not guaranteed to match the node ordering
# during the run, which itself is not significant. The columns contain field
# values in whatever units where used in the simulation. Usually MKS.
# Columns are:
#  | id | eos_id | x | y | z | vx | vy | vz | m | rho | p | T | U | hmin | hmax |
#
# Column legend:
#    
#        id - an integer identifier of the node list this node came from
#    eos_id - an integer identifier of the material eos used with this node list
#     x,y,z - node space coordinates 
#  vx,vy,vz - node velocity components
#         m - node mass
#       rho - mass density
#         p - pressure
#         T - temperature
#         U - specific internal energy
# hmin,hmax - smallest and largest half-axes of the smoothing ellipsoid 
#
# Tip: load table into python with np.loadtxt()
#
################################################################################
"""

global material_dictionary
# A dictionary of unique short tags for commonly used material EOSs.
# We use this in spite of the added complexity to allow users of uss to specify 
# nothing more than a unique string material "tag" as a complete choice of eos.
# All uss runs should use this tag and the construct_eos_for_material method 
# instead of calling the spheral eos constructors directly. This also allows 
# keeping a record of what material was used in each run, and thus allows a hands
# free importing of precooked planets into new runs.
# IMPORTANT: add new entries AFTER old ones to preserve unique ids.
material_dictionary = {}

material_dictionary['water'] = dict(
        eos_type = 'Tillotson',
        eos_constructor = sph.TillotsonEquationOfState,
        eos_arguments = {'materialName':'water'},
        eos_id = len(material_dictionary.keys()) + 1,
        )

material_dictionary['h2oice'] = dict(
        eos_type = 'Tillotson',
        eos_constructor = sph.TillotsonEquationOfState,
        eos_arguments = {'materialName':'pure ice'},
        eos_id = len(material_dictionary.keys()) + 1,
        )

material_dictionary['dirtyice'] = dict(
        eos_type = 'Tillotson',
        eos_constructor = sph.TillotsonEquationOfState,
        eos_arguments = {'materialName':'30% silicate ice'},
        eos_id = len(material_dictionary.keys()) + 1,
        )

material_dictionary['granite'] = dict(
        eos_type = 'Tillotson',
        eos_constructor = sph.TillotsonEquationOfState,
        eos_arguments = {'materialName':'granite'},
        eos_id = len(material_dictionary.keys()) + 1,
        )

material_dictionary['basalt'] = dict(
        eos_type = 'Tillotson',
        eos_constructor = sph.TillotsonEquationOfState,
        eos_arguments = {'materialName':'basalt'},
        eos_id = len(material_dictionary.keys()) + 1,
        )

material_dictionary['nylon'] = dict(
        eos_type = 'Tillotson',
        eos_constructor = sph.TillotsonEquationOfState,
        eos_arguments = {'materialName':'nylon'},
        eos_id = len(material_dictionary.keys()) + 1,
        )

material_dictionary['sio2'] = dict(
        eos_type = 'M/ANEOS',
        eos_constructor = None,
        eos_arguments = {},
        eos_id = len(material_dictionary.keys()) + 1,
        )
# End material_dictionary
