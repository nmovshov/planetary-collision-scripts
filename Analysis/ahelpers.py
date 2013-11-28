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

def _test():
    print "alo"
    fnl=FNLMeta()
    print fnl.nb_columns
    pass

