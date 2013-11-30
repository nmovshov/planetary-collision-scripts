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
import matplotlib as plt

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
    """Load node list data from file and parse out to a struct."""

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

    fnl = FNLData()
    fnl.id =   data[:,  FNLMeta.nl_id_col]
    fnl.eos =  data[:, FNLMeta.eos_id_col]
    fnl.x =    data[:,      FNLMeta.x_col]
    fnl.y =    data[:,      FNLMeta.y_col]
    fnl.z =    data[:,      FNLMeta.z_col]
    fnl.vx =   data[:,     FNLMeta.vx_col]
    fnl.vy =   data[:,     FNLMeta.vy_col]
    fnl.vz =   data[:,     FNLMeta.vz_col]
    fnl.m =    data[:,      FNLMeta.m_col]
    fnl.rho =  data[:,    FNLMeta.rho_col]
    fnl.P =    data[:,      FNLMeta.P_col]
    fnl.T =    data[:,      FNLMeta.T_col]
    fnl.U =    data[:,      FNLMeta.U_col]
    fnl.hmin = data[:,   FNLMeta.hmin_col]
    fnl.hmax = data[:,   FNLMeta.hmax_col]

    return fnl

def _test():
    print "alo"
    pass

