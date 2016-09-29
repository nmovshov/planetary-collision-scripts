#-------------------------------------------------------------------------------
#   Analysis Helpers - A collection of some convenience functions for reuse in
#                      automated or interactive analysis of collision runs.
#
# Author: nmovshov at gmail dot com
#-------------------------------------------------------------------------------
import sys, os, shutil
import numpy as np

class FNLMeta:
    """A struct with info about the layout of .fnl files."""

    nb_columns = 15
    nl_id_col = 0
    eos_id_col = 1
    x_col = 2
    y_col = 3
    z_col = 4
    vx_col = 5
    vy_col = 6
    vz_col = 7
    m_col = 8
    rho_col = 9
    P_col = 10
    T_col = 11
    U_col = 12
    hmin_col = 13
    hmax_col = 14
    pass

class FNLData:
    """An empty struct that can be used to hold essential node list data."""
    pass

def load_fnl(filename):
    """Load node list data from file and parse out to a struct.
    
    The file filename is assumed to contain data from one or more node lists that
    have been flattened using shelpers.pflatten_node_list_list. Minimal checking
    is applied but I assume a responsible user. A list of flattened node lists
    will have one or more node lists identified by consecutive, integer, zero-
    based id. This method will return a single FNLData struct with combining all
    records in the file.
    """

    # Read raw data
    assert isinstance(filename, str)
    try:
        data = np.loadtxt(filename)
    except:
        raise StandardError(
            "ERROR: Could not read data from file {}".format(filename))
    if (data.ndim != 2) or (data.shape[1] != FNLMeta.nb_columns):
        raise StandardError(
            "ERROR: {} does not appear to be a valid flattened node list.".format(
            filename))

    fnl = FNLData()

    fnl.id   = data[:,  FNLMeta.nl_id_col]
    fnl.eos  = data[:, FNLMeta.eos_id_col]
    fnl.x    = data[:,      FNLMeta.x_col]
    fnl.y    = data[:,      FNLMeta.y_col]
    fnl.z    = data[:,      FNLMeta.z_col]
    fnl.vx   = data[:,     FNLMeta.vx_col]
    fnl.vy   = data[:,     FNLMeta.vy_col]
    fnl.vz   = data[:,     FNLMeta.vz_col]
    fnl.m    = data[:,      FNLMeta.m_col]
    fnl.rho  = data[:,    FNLMeta.rho_col]
    fnl.P    = data[:,      FNLMeta.P_col]
    fnl.T    = data[:,      FNLMeta.T_col]
    fnl.U    = data[:,      FNLMeta.U_col]
    fnl.hmin = data[:,   FNLMeta.hmin_col]
    fnl.hmax = data[:,   FNLMeta.hmax_col]
    fnl.nbNodes = len(data)
    fnl.r = np.hypot(fnl.x, np.hypot(fnl.y, fnl.z))

    return fnl

def pack_fnl(data):
    """Pack fnl array to fnl struct."""

    # Minimal input control
    assert isinstance(data, np.ndarray)
    assert data.ndim == 2 and data.shape[1] == FNLMeta.nb_columns

    # Pack fnl
    fnl = FNLData()
    fnl.id   = data[:,  FNLMeta.nl_id_col]
    fnl.eos  = data[:, FNLMeta.eos_id_col]
    fnl.x    = data[:,      FNLMeta.x_col]
    fnl.y    = data[:,      FNLMeta.y_col]
    fnl.z    = data[:,      FNLMeta.z_col]
    fnl.vx   = data[:,     FNLMeta.vx_col]
    fnl.vy   = data[:,     FNLMeta.vy_col]
    fnl.vz   = data[:,     FNLMeta.vz_col]
    fnl.m    = data[:,      FNLMeta.m_col]
    fnl.rho  = data[:,    FNLMeta.rho_col]
    fnl.P    = data[:,      FNLMeta.P_col]
    fnl.T    = data[:,      FNLMeta.T_col]
    fnl.U    = data[:,      FNLMeta.U_col]
    fnl.hmin = data[:,   FNLMeta.hmin_col]
    fnl.hmax = data[:,   FNLMeta.hmax_col]
    fnl.nbNodes = len(data)
    fnl.r = np.hypot(fnl.x, np.hypot(fnl.y, fnl.z))

    # Return
    return fnl

def unpack_fnl(fnl):
    """Unpack fnl struct to rectangular array."""

    # Minimal input control
    assert isinstance(fnl, FNLData)

    # Allocate array
    data = np.zeros([fnl.nbNodes,FNLMeta.nb_columns])

    # Unpack fnl struct
    data[:,  FNLMeta.nl_id_col] = fnl.id
    data[:, FNLMeta.eos_id_col] = fnl.eos
    data[:,      FNLMeta.x_col] = fnl.x
    data[:,      FNLMeta.y_col] = fnl.y
    data[:,      FNLMeta.z_col] = fnl.z
    data[:,     FNLMeta.vx_col] = fnl.vx
    data[:,     FNLMeta.vy_col] = fnl.vy
    data[:,     FNLMeta.vz_col] = fnl.vz
    data[:,      FNLMeta.m_col] = fnl.m
    data[:,    FNLMeta.rho_col] = fnl.rho
    data[:,      FNLMeta.P_col] = fnl.P
    data[:,      FNLMeta.T_col] = fnl.T
    data[:,      FNLMeta.U_col] = fnl.U
    data[:,   FNLMeta.hmin_col] = fnl.hmin
    data[:,   FNLMeta.hmax_col] = fnl.hmax

    # Return
    return data

def save_fnl(filename, fnl, header=None):
    """Unpack fnl struct to array and save to ascii file."""

    # Minimal input control
    assert isinstance(fnl, FNLData)
    assert isinstance(filename, str)

    # Allocate array
    data = np.zeros([fnl.nbNodes,FNLMeta.nb_columns])

    # Unpack fnl struct
    data[:,  FNLMeta.nl_id_col] = fnl.id
    data[:, FNLMeta.eos_id_col] = fnl.eos
    data[:,      FNLMeta.x_col] = fnl.x
    data[:,      FNLMeta.y_col] = fnl.y
    data[:,      FNLMeta.z_col] = fnl.z
    data[:,     FNLMeta.vx_col] = fnl.vx
    data[:,     FNLMeta.vy_col] = fnl.vy
    data[:,     FNLMeta.vz_col] = fnl.vz
    data[:,      FNLMeta.m_col] = fnl.m
    data[:,    FNLMeta.rho_col] = fnl.rho
    data[:,      FNLMeta.P_col] = fnl.P
    data[:,      FNLMeta.T_col] = fnl.T
    data[:,      FNLMeta.U_col] = fnl.U
    data[:,   FNLMeta.hmin_col] = fnl.hmin
    data[:,   FNLMeta.hmax_col] = fnl.hmax

    # Write to file
    if header is None:
        header = fnl_header_default.format(fnl.nbNodes)
    else:
        assert isinstance(header, str)
    format = 2*['%2d'] + (FNLMeta.nb_columns - 2)*['%12.5e']
    np.savetxt(filename, data, header=header, fmt=format)

    # Return
    return

def load_multi_fnl(filename):

    """Load node list data from file and parse out to a struct.
    
    The file filename is assumed to contain data from one or more node lists that
    have been flattened using shelpers.pflatten_node_list_list. Minimal checking
    is applied but I assume a responsible user. A list of flattened node lists
    will have one or more node lists identified by consecutive, integer, zero-
    based id. This method will return a tuple of FNLData structs with convenient
    field names. If you know in advance the number of node lists in the file you
    can use tuple unpacking to get individual FNLData structs. In the case of a
    file containing a single node list a single struct is returned instead of a
    tuple.
    """

    # Read raw data
    assert isinstance(filename, str)
    try:
        data = np.loadtxt(filename)
    except:
        raise StandardError("ERROR: Could not read data from file {}".format(
            filename))
    if (data.ndim != 2) or (data.shape[1] != FNLMeta.nb_columns):
        raise StandardError(
            "ERROR: {} does not appear to be a valid flattened node list.".format(
            filename))

    nbLists = int(max(data[:,0])) + 1
    fnl = tuple([FNLData() for k in range(nbLists)])

    for k in range(nbLists):
        kmask = data[:,FNLMeta.nl_id_col]==k

        fnl[k].id =   data[kmask,  FNLMeta.nl_id_col]
        fnl[k].eos =  data[kmask, FNLMeta.eos_id_col]
        fnl[k].x =    data[kmask,      FNLMeta.x_col]
        fnl[k].y =    data[kmask,      FNLMeta.y_col]
        fnl[k].z =    data[kmask,      FNLMeta.z_col]
        fnl[k].vx =   data[kmask,     FNLMeta.vx_col]
        fnl[k].vy =   data[kmask,     FNLMeta.vy_col]
        fnl[k].vz =   data[kmask,     FNLMeta.vz_col]
        fnl[k].m =    data[kmask,      FNLMeta.m_col]
        fnl[k].rho =  data[kmask,    FNLMeta.rho_col]
        fnl[k].P =    data[kmask,      FNLMeta.P_col]
        fnl[k].T =    data[kmask,      FNLMeta.T_col]
        fnl[k].U =    data[kmask,      FNLMeta.U_col]
        fnl[k].hmin = data[kmask,   FNLMeta.hmin_col]
        fnl[k].hmax = data[kmask,   FNLMeta.hmax_col]
        fnl[k].nbNodes = sum(kmask)
        fnl[k].r = np.hypot(fnl[k].x,np.hypot(fnl[k].y,fnl[k].z))

    if len(fnl)>1:
        return fnl
    else:
        return fnl[0]

def plot_P_vs_r(fnl, bblock=False):
    """Plot pressure of nodes against distance from origin."""

    import matplotlib as mpl
    import matplotlib.pyplot as plt

    assert isinstance(fnl,(FNLData,tuple))
    if isinstance(fnl,FNLData):
        fnl = (fnl,)

    fig = plt.figure()
    axe = plt.axes()
    plt.xlabel('Radius [km]')
    plt.ylabel('Pressure [GPa]')
    plt.grid()
    for nl in fnl:
        assert isinstance(nl,FNLData)
        x = np.sort(nl.r)
        y = nl.P[np.argsort(nl.r)]
        plt.plot(x/1e3, y/1e9)
        pass
    plt.show(block=bblock)
    return (fig,axe)

def plot_rho_vs_r(fnl, bblock=False):
    """Plot density of nodes against distance from origin."""

    import matplotlib as mpl
    import matplotlib.pyplot as plt

    assert isinstance(fnl,(FNLData,tuple))
    if isinstance(fnl,FNLData):
        fnl = (fnl,)

    fig = plt.figure()
    axe = plt.axes()
    plt.xlabel('Radius [km]')
    plt.ylabel('Mass density [kg/m^3]')
    plt.grid()
    for nl in fnl:
        assert isinstance(nl,FNLData)
        x = np.sort(nl.r)
        y = nl.rho[np.argsort(nl.r)]
        plt.plot(x/1e3, y)
        pass
    plt.show(block=bblock)
    return (fig,axe)

def plot_P_vs_r_output(dirname='.', bblock=False):
    """Plot P(r) for all fnl files in a directory."""

    import matplotlib as mpl
    import matplotlib.pyplot as plt

    assert isinstance(dirname,str)
    assert os.path.isdir(dirname)

    fnl_files = []
    for root, dirs, files in os.walk(dirname):
        fnl_files += [os.path.join(root,fn) for fn in files if 
                                                fn.endswith(('.fnl','.fnl.gz'))]
    fnl_files.sort()
    if len(fnl_files) == 0:
        print "No .fnl or .fnl.gz files found in directory."
        return
    all_fnls = [load_multi_fnl(f) for f in fnl_files]

    fig = plt.figure()
    nb_rows = np.ceil(np.sqrt(len(all_fnls)))
    nb_cols = np.ceil(len(all_fnls)/nb_rows)
    for k in range(len(all_fnls)):
        fnl = all_fnls[k]
        if isinstance(fnl,FNLData):
            fnl = (fnl,)
        plt.subplot(nb_rows,nb_cols,k+1)
        for nl in fnl:
            assert isinstance(nl,FNLData)
            x = np.sort(nl.r)
            y = nl.P[np.argsort(nl.r)]
            plt.plot(x/1e3, y/1e9)
            pass
        pass

    plt.show(block=bblock)    
    return fig

def plot_rho_vs_r_output(dirname='.', bblock=False):
    """Plot rho(r) for all fnl files in a directory."""

    import matplotlib as mpl
    import matplotlib.pyplot as plt

    assert isinstance(dirname,str)
    assert os.path.isdir(dirname)

    fnl_files = []
    for root, dirs, files in os.walk(dirname):
        fnl_files += [os.path.join(root,fn) for fn in files if 
                                                fn.endswith(('.fnl','.fnl.gz'))]
    fnl_files.sort()
    if len(fnl_files) == 0:
        print "No .fnl or .fnl.gz files found in directory."
        return
    all_fnls = [load_multi_fnl(f) for f in fnl_files]

    fig = plt.figure()
    nb_rows = np.ceil(np.sqrt(len(all_fnls)))
    nb_cols = np.ceil(len(all_fnls)/nb_rows)
    for k in range(len(all_fnls)):
        fnl = all_fnls[k]
        if isinstance(fnl,FNLData):
            fnl = (fnl,)
        plt.subplot(nb_rows,nb_cols,k+1)
        for nl in fnl:
            assert isinstance(nl,FNLData)
            x = np.sort(nl.r)
            y = nl.rho[np.argsort(nl.r)]
            plt.plot(x/1e3, y)
            pass
        pass

    plt.show(block=bblock)    
    return fig

def plot_XY_scatter(fnl, bblock=False):
    """XY scatter plot of node positions with color density."""

    import matplotlib as mpl
    import matplotlib.pyplot as plt

    assert isinstance(fnl,(FNLData,tuple))
    if isinstance(fnl,FNLData):
        fnl = (fnl,)

    fig = plt.figure()
    axe = plt.axes()
    plt.xlabel('X [km]')
    plt.ylabel('Y [km]')
    plt.grid()
    for nl in fnl:
        assert isinstance(nl,FNLData)
        x = nl.x/1e3
        y = nl.y/1e3
        rho = nl.rho
        h = nl.hmin
        plt.scatter(x, y, s=20, c=rho)
        pass
    plt.show(block=bblock)
    return (fig,axe)

def ejectify_fnl(fnl, method='naor1'):
    """Quick-and-dirty detection of ejecta field from SPHERAL output fnl."""

    # Minimal input control
    assert isinstance(fnl, FNLData)
    assert method in ('naor1', 'jutzi')

    # Extract values from loaded fnl
    data = unpack_fnl(fnl)
    pos = data[:,FNLMeta.x_col:FNLMeta.z_col+1]
    vel = data[:,FNLMeta.vx_col:FNLMeta.vz_col+1]
    m = data[:,FNLMeta.m_col]

    # Detect largest bound mass to serve as primary
    import bound_mass as bm
    [M, ind] = bm.bound_mass(pos, vel, m, method=method)

    # Move to primary rest CoM frame
    X = [sum(pos[ind,0]*m[ind]), sum(pos[ind,1]*m[ind]), sum(pos[ind,2]*m[ind])]/M
    V = [sum(vel[ind,0]*m[ind]), sum(vel[ind,1]*m[ind]), sum(vel[ind,2]*m[ind])]/M
    data[:,FNLMeta.x_col:FNLMeta.z_col+1] = pos - X
    data[:,FNLMeta.vx_col:FNLMeta.vz_col+1] = vel - V

    # Remove primary and repack to fnl
    data = data[~ind]
    ejecta = pack_fnl(data)

    return ejecta, ~ind

def _test():
    print "alo"
    pass

global fnl_header_default
fnl_header_default = """\
###############################################################################
 This file contains output data from a Spheral++ simulation, including all
 field variables as well as some diagnostic data and node meta data. This
 file should contain {} data lines, one per SPH node used in the simulation.
 Line order is not significant and is not guaranteed to match the node ordering
 during the run, which itself is not significant. The columns contain field
 values in whatever units where used in the simulation. Usually MKS.
 Columns are:
  | id | eos_id | x | y | z | vx | vy | vz | m | rho | p | T | U | hmin | hmax |

 Column legend:

        id - an integer identifier of the node list this node came from
    eos_id - an integer identifier of the material eos used with this node list
     x,y,z - node space coordinates
  vx,vy,vz - node velocity components
         m - node mass
       rho - mass density
         p - pressure
         T - temperature
         U - specific internal energy
 hmin,hmax - smallest and largest half-axes of the smoothing ellipsoid

 Tip: load table into python with np.loadtxt()

###############################################################################
"""

global fnl_header_ejecta
fnl_header_ejecta = """\
###############################################################################
 This file contains the ejecta field from an impact simulation. Included are
 field variables as well as some diagnostic data and node meta data. Original
 simulation had M={:0.4g} kg in {} nodes.
 This file was generated by identifying the largest gravitationally bound
 fragment, with M_lb={:0.4g} kg, moving to that fragment's center of mass
 frame, and then removing all nodes belonging to the largest fragment.
 This file should contain {} data lines, one per SPH node identified
 as ejecta. Line order is not significant. The columns contain field values in
 whatever units where used in the simulation, hopefully MKS.
 Columns are:
 | id | eos_id | x | y | z | vx | vy | vz | m | rho | p | T | U | h | R_ini |

 Column legend:

        id - an integer identifier of the node list this node came from
    eos_id - an integer identifier of the material eos used with this node list
     x,y,z - node space coordinates
  vx,vy,vz - node velocity components
         m - node mass
       rho - mass density
         p - pressure
         T - temperature
         U - specific internal energy
         h - smoothing length
     R_ini - pre-impact radius from center of parent body (normalized)

 Tip: load table into python with np.loadtxt()
 Note on eos_id: 2=h2oice (tillotson), 5=basalt (tillotson)

###############################################################################
"""

if __name__ == "__main__":
    _test()
    pass
