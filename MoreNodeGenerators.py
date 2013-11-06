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
import numpy as np
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
        assert nLayers > 1
        assert type(rho)==type(rMin)==type(rMax)==type(nNodePerh)==float
        assert rho > 0.0
        assert rMax > rMin >= 0.0
        assert nNodePerh >= 1.0

        # Store key parameters in the generator object.
        self.nLayers = nLayers
        self.rho = rho # we'll convert this to a list later.
        self.rMin = rMin
        self.rMax = rMax
        self.nNodePerh = nNodePerh

        # Create the required lists (empty).
        self.x=[]
        self.y=[]
        self.z=[]
        self.m=[]
        self.H=[]

        # Fill lists with calculate positions, masses, Hs.
        self._generate_equally_spaced_shells()

        # Make rho into a list.
        self.rho = [self.rho] * len(self.m)

        # Have the base class break up the serial node distribution
        # for parallel cases.
        NodeGeneratorBase.__init__(self, True,
                                   self.x, self.y, self.z, self.m, self.H)
        pass

    #---------------------------------------------------------------------------
    # The actual generator algorithm
    #---------------------------------------------------------------------------
    def _generate_equally_spaced_shells(self):
        """Given shell spacing, fill sphere with equally spaced nodes.

        The idea is simple. The requested number of equally spaced shells defines
        the linear spacing between nodes, dl. Use that linear spacing to fill the
        sphere by filling up slices first, then stacks, then shells.
        
        The nodes are places in the center of the shells, so that if rMin and
        rMax specify the boundaries of the sphere, the node at rMin is the center
        of the shell extending from rMin-dl/2 to rMin+dl/2, and similarly for
        rMax. (Therefor the "surface" of the planet would be at rMax+dl/2.)

        The stacks are built by placing a node at the "north" and "south" poles,
        and as many equally spaced nodes as needed in between such that the
        linear distance between nodes on the great circle is close to dl.

        Then each stack is divided to equally spaced slices. But the slices don't
        all being at a "prime meridian." Instead each is shifted a linear dl from
        the last.

        TODO: The volume element associated with each node should be calculated.
        at the moment, all nodes are given an equal mass, but this is a rough
        approximation.
        """
        dl = (self.rMax-self.rMin)/(self.nLayers-1) # Constant linear spacing.

        # Fill up shells...
        shells = np.linspace(self.rMin,self.rMax,self.nLayers)
        for r in shells:
            if r==0.0:
                self.x.append(0.0)
                self.y.append(0.0)
                self.z.append(0.0)
                continue
            # With stacks...
            dG = dl/r
            nGs = int(pi/dG)+1 # yes, +1, linspace includes end points
            stacks = np.linspace(0.0,pi,nGs)
            for G in stacks:
                if G==0.0 or G==pi:
                    self.x.append(0.0)
                    self.y.append(0.0)
                    self.z.append(r*cos(G))
                    continue
                # Made of slices.
                dq = dl/(r*sin(G))
                nqs = int(2*pi/dq)+1 # yes, +1, linspace includes end points
                slices = np.linspace(0.0,2*pi,nqs)[1:] # 0=2pi
                slices += np.random.uniform(0,pi/4)
                for q in slices:
                    self.x.append(r*sin(G)*cos(q))
                    self.y.append(r*sin(G)*sin(q))
                    self.z.append(r*cos(G))
                    pass
                pass
            pass

        # Assign smoothing tensor and mass.
        h0 = 1.0/(dl*self.nNodePerh)
        nominalH = SymTensor3d(h0, 0.0, 0.0,
                               0.0, h0, 0.0,
                               0.0, 0.0, h0)
        self.H = [nominalH]*len(self.x)
        
        #TODO: improve this!
        totalMass = self.rho*4.0*pi/3.0 * (self.rMax**3 - self.rMin**3)
        nominalMassPerNode = totalMass / len(self.x)
        self.m = [nominalMassPerNode]*len(self.x)
        
        # And Bob's our uncle.
        pass # end of method.

    #---------------------------------------------------------------------------
    # Required methods from NodeGeneratorBase
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        assert i >= 0 and i < len(self.x)
        assert len(self.x) == len(self.y) == len(self.z)
        return Vector3d(self.x[i], self.y[i], self.z[i])

    def localMass(self, i):
        assert i >= 0 and i < len(self.m)
        return self.m[i]

    def localMassDensity(self, i):
        assert i >= 0 and i < len(self.rho)
        return self.rho[i]

    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]
