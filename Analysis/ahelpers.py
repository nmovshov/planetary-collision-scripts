#-------------------------------------------------------------------------------
#   Analysis Helpers - A collection of some convenience functions for reuse in
#                      automated or interactive analysis of collision runs.
#
# Author: nmovshov at gmail dot com
#-------------------------------------------------------------------------------
import sys, os, shutil
import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt

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
    is applied, but I assume a responsible user. A list of flattened node lists
    will have one or more node lists identified by consecutive, integer, zero-
    based id. This method will return a tuple of FNLData structs, with convenient
    field names. If you know in advance the number of node lists in the file, you
    can use tuple unpacking to get individual FNLData structs. In the case of a
    file containing a single node list, a single struct is returned instead of a
    tuple.
    """

    # Read raw data
    assert isinstance(filename, str)
    try:
        data = np.loadtxt(filename)
    except:
        print "ERROR: Could not read data from file {}".format(filename)
        return None
    if (data.ndim != 2) or (data.shape[1] != FNLMeta.nb_columns):
        print "ERROR: {} does not appear to be a valid flattened node list.".format(
                filename)
        return None

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

def plot_P_vs_r(fnl):
    """Plot pressure of nodes against distance from origin."""

    assert isinstance(fnl,(FNLData,tuple))
    if isinstance(fnl,FNLData):
        fnl = (fnl,)

    fig = plt.figure()
    axe = plt.axes()
    plt.xlabel('Radius [m]')
    plt.ylabel('Pressure [GPa]')
    plt.grid()
    for nl in fnl:
        assert isinstance(nl,FNLData)
        plt.plot(nl.r,nl.P/1e9,'.')
        pass
    plt.show(block=False)
    return (fig,axe)

def plot_P_vs_r_output(dirname):
    """Plot P(r) for all fnl files in a directory."""

    assert isinstance(dirname,str)
    assert os.path.isdir(dirname)

    fnl_files = []
    for root, dirs, files in os.walk(dirname):
        fnl_files += [os.path.join(root,fn) for fn in files if 
                                                fn.endswith(('.fnl','.fnl.gz'))]
    fnl_files.sort()
    all_fnls = [load_fnl(f) for f in fnl_files]

    fig = plt.figure()
    nb_rows = np.ceil(np.sqrt(len(all_fnls)))
    nb_cols = np.ceil(len(all_fnls)/nb_rows)
    for k in range(len(all_fnls)):
        fnl = all_fnls[k]
        if isinstance(fnl,FNLData):
            fnl = (fnl,)
        plt.subplot(nb_rows,nb_cols,k+1)
        for nl in fnl:
            plt.plot(nl.r,nl.P/1e9,'.')
            pass
        pass

    plt.show(block=False)    
    return fig

def _test():
    print "alo"
    pass

