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
                 nNodePerh = 2.01,
                 EOS = None,):

        # Some assertions for convenience. Not supposed to be an airtight seal.
        assert type(nLayers) == int
        assert nLayers > 1
        assert type(rho)==type(rMin)==type(rMax)==type(nNodePerh) == float
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
        
        # I use a volume field to facilitate modifying mass post construction.
        self.V=[] # Volume associated with node.

        # Fill lists with calculated positions, masses, Hs.
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
        linear distance between nodes on the great circle is close to dl. These
        colatitudes mark the _center_ of the stack.

        Each stack is then divided to equally spaced slices. But the slices don't
        all begin at a "prime meridian." Instead each is shifted a linear dl from
        the last.

        TODO: The volume element associated with each node should be calculated.
        at the moment, all nodes are given an equal mass, but this is a rough
        approximation. The V field is not yet implemented.
        """
        dl = (self.rMax-self.rMin)/(self.nLayers-1) # Constant linear spacing.

        # Fill up shells...
        dr = dl
        shells = np.linspace(self.rMin,self.rMax,self.nLayers)
        for r in shells:
            if r==0.0:
                self.x.append(0.0)
                self.y.append(0.0)
                self.z.append(0.0)
                continue
            # With stacks...
            dG = dl/r
            nGs = int(pi/dG)
            stacks = np.linspace(0.0,pi,nGs)
            for G in stacks:
                if G==0.0 or G==pi:
                    self.x.append(0.0)
                    self.y.append(0.0)
                    self.z.append(r*cos(G))
                    continue
                # Made of slices.
                dq = dl/(r*sin(G))
                nqs = int(2*pi/dq)
                slices = np.linspace(0.0,2*pi,nqs+1)[1:] # 0=2pi
                #slices += np.random.uniform(0,pi/4)
                G0=G-dG/2.0
                G1=G+dG/2.0
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


class HexagonalClosePacking(NodeGeneratorBase):
    """Put nodes on HCP lattice and optionally chisel out a sphere."""

    #---------------------------------------------------------------------------
    # The constructor
    #---------------------------------------------------------------------------
    def __init__(self, nx, rho,
                 xMin = 0.0,
                 xMax = 1.0,
                 rMin = 0.0,
                 rMax = 1e200,
                 nNodePerh = 2.01,
                 EOS = None,):
        """Class constructor for the HCP lattice node generator. Note that the
           generated coordinates are those of the CENTER of the node. Also note
           that because the lines and layers of the HCP lattice do not end at
           equal coordinate limits, the volume enclosed by the lattice is only
           approximately equal to (xMax-xMin)**3. 
           
          Parameters
          ----------
          nx : int > 0
              Number of nodes across domain. There are no corresponding ny or nz
              because we want to ensure a cubic lattice, with equally spaced nodes
              in all directions.
          rho : float > 0
              Density used to assign node masses. For now, this is a constant.
          xMin : float, optional
              Left edge of lattice. Default is 0.0.
          xMax : float, optional
              Right edge of lattice. Default is 1.0.
          rMin : float >=0, optional
              After lattice is built, nodes whose distance from the center of the
              lattice is less than rMin will be culled. Default is 0.0.
          rMax : float > 0, optional
              After lattice is built, nodes whose distance from the center of the
              lattice is greater than rMax will be culled. Default is 1e200.
          nNodePerh : float > 1.0, optional
              Nodes are assigned an inverse smoothing length of 1/(d*nNodePerh),
              where d is the lattice spacing. Default is 2.01.
          EOS : Spheral.SolidMaterial.SolidEquationOfState, optional.
              Place holder, for future use. Default is None.
        """

        # Some assertions for convenience. Not supposed to be an airtight seal.
        assert isinstance(nx,int) and nx > 1
        assert isinstance(rho,float) and rho > 0.0
        assert isinstance(xMin,float) and isinstance(xMax,float)
        assert xMax > xMin
        assert isinstance(rMin,float) and isinstance(rMax,float)
        assert rMax > rMin >= 0.0
        assert isinstance(nNodePerh,float) and nNodePerh >= 1.0
        assert True # place holder for eos

        # Store key parameters in the generator object.
        self.nx = nx
        self.rho = rho
        self.xMin = xMin
        self.xMax = xMax
        self.rMin = rMin
        self.rMax = rMax
        self.nNodePerh = nNodePerh
        self.EOS = EOS

        # Create the required lists (empty).
        self.x=[]
        self.y=[]
        self.z=[]
        self.m=[]
        self.H=[]

        # I use a volume field to facilitate modifying mass post construction.
        self.lattice_spacing = (self.xMax-self.xMin)/(self.nx)
        self.lattice_volume = (self.xMax-self.xMin)**3
        self.V = []

        # Fill lists with calculated positions, masses, Hs.
        self._generate_hcp_lattice()

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
    def _generate_hcp_lattice(self):
        """Here we generate hcp lattice positions. Each node will have 12 nearest
           neighbors, all a distance d from it, where d, the lattice spacing, is
           the diameter of the spheres that can be placed on the lattice in
           closest packing.
        """

        d = self.lattice_spacing
        r = d/2
        xstep = d
        ystep = sqrt(3)/2 * d
        zstep = sqrt(6)/3 * d
        nx = self.nx
        ny = int(nx * xstep/ystep) + 1
        nz = int(nx * xstep/zstep) + 1
        
        for k in range(nz):
            z = r + k*zstep
            for j in range(ny):
                if np.mod(k,2)==0:
                    y = r + j*ystep
                else:
                    y = r + d*sqrt(3)/6 + j*ystep
                for i in range(nx):
                    if np.mod(j,2)==0:
                        x = r + r*np.mod(k,2) + i*xstep
                    else:
                        x = d - r*np.mod(k,2) + i*xstep
                    self.x.append(x)
                    self.y.append(y)
                    self.z.append(z)

                    pass
                pass
            pass

        # Assign masses and smoothing tensors.
        nominalCellVolume = self.lattice_volume/len(self.x)
        self.V = [nominalCellVolume] * len(self.x)

        nominalCellMass = nominalCellVolume*self.rho
        self.m = [nominalCellMass] * len(self.x)

        h0 = 1.0/(d*self.nNodePerh)
        nominalH = SymTensor3d(h0,  0.0, 0.0,
                               0.0, h0,  0.0,
                               0.0, 0.0, h0)
        self.H = [nominalH] * len(self.x)

        # And Bob's our uncle.
        pass

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
