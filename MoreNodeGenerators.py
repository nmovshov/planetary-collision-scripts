#-------------------------------------------------------------------------------
# A few more node generators, dealing with some special cases.
#-------------------------------------------------------------------------------
from math import *
from NodeGeneratorBase import *
from Spheral import Vector3d
from Spheral import Tensor3d
from Spheral import SymTensor3d
from Spheral import vector_of_int, vector_of_double 
from Spheral import vector_of_SymTensor3d, vector_of_vector_of_double
import mpi
procID = mpi.rank
nProcs = mpi.procs

class EqualSpacingSphericalShells(NodeGeneratorBase):
    """Aiming for spherically symmetric distribution of equally spaced nodes."""

    #---------------------------------------------------------------------------
    # The constructor
    #---------------------------------------------------------------------------
    def __init__(self, nLayers, rho,
                 rMin = 0.0,
                 rMax = 1.0,
                 nNodePerh = 2.01,):

        # Some assertions for convenience. Not supposed to be an airtight seal.
        assert type(nLayers) == int
        assert nLayers > 0
        assert type(rho)==type(rMin)==type(rMax)==type(nNodePerh)==float
        assert rho > 0.0
        assert rMax > rMin >= 0.0
        assert nNodePerh >= 1.0

        # Store key parameters in the generator object.
        self.nLayers = nLayers
        self.rho = rho # we'll convert this to a list later
        self.rMin = rMin
        self.rMax = rMax
        self.nNodePerh = nNodePerh

        # Generate lists for position, mass, and H
        (self.x, self.y, self.z, self.m, self.H) = \
                self.generate_equally_spaced_shells()

        # Make rho into a list.
        self.rho = [self.rho] * len(self.m)

        # Have the base class break up the serial node distribution
        # for parallel cases.
        NodeGeneratorBase.__init__(self, True,
                                   self.x, self.y, self.z, self.m, self.H)
        print "alo" 
        pass

    #---------------------------------------------------------------------------
    # The actual generator algorithm
    #---------------------------------------------------------------------------
    def generate_equally_spaced_shells(self):
        x,y,z,m,H=[1],[2],[3],[4],[5]
        return x, y, z, m, H
        pass

    #---------------------------------------------------------------------------
    # Required methods from NodeGeneratorBase
    #---------------------------------------------------------------------------
