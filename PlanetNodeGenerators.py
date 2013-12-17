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
    """Put nodes on spherical shells keeping equal nearest-neighbor spacing."""

    #---------------------------------------------------------------------------
    # The constructor
    #---------------------------------------------------------------------------
    def __init__(self, nLayers, rho,
                 rMin = 0.0,
                 rMax = 1.0,
                 nNodePerh = 2.01,
                 EOS = None,):
        """Class constructor for the spherical shells node generator. 

        The generated coordinates are those of the CENTER of the node. So, for
        example, no node will be at rMin or rMax. Also note that the volume of
        space associated with each node is NOT equal. For low resolution runs,
        this can result in large mass differences between nodes. The effect of
        that is unknown though.

        Parameters
        ----------
        nLayers : int > 1
            Number of layers, or shells. The use does not specify a number of
            stacks or slices because these will be determined by the equal
            spacing requirement. Note: while this parameter is basically the 
            one controlling run resolution, its meaning is different than the
            nx parameter for cubic grids, and the resulting number of nodes
            is usually higher.
        rho : float > 0
            Density used to assign node masses. For now, this is a constant.
        rMin : float >=0, optional
            Inner edge of domain. NOT the coordinate of innermost node!
            Default is 0.0.
        rMax : float > 0, optional
            Outermost edge of domain. NOT the coordinate of outermost node!
            Default is 1.0.
        nNodePerh : float > 1.0, optional
            Nodes are assigned an inverse smoothing length of 1/(dl*nNodePerh),
            where dl is the calculated node spacing. Default is 2.01.
        EOS : Spheral.SolidMaterial.SolidEquationOfState, optional.
            Place holder, for future use. Default is None.
        """

        # Some assertions for convenience. Not supposed to be an airtight seal.
        assert isinstance(nLayers,int) and nLayers > 1
        assert isinstance(rho,float) and rho > 0.0
        assert isinstance(rMin,float) and isinstance(rMax,float)
        assert rMax > rMin >= 0.0
        assert isinstance(nNodePerh,float) and nNodePerh >= 1.0

        # Store key parameters in the generator object.
        self.nLayers = nLayers
        self.rho = rho # will become a list when we get the node count.
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
        
        # Convenience and diagnostic fields. (Don't overdo.)
        self.linear_spacing = 0.0
        self.domain_volume = 4*pi/3*rho*(rMax**3-rMin**3)
        self.V = [] 

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
        """Fill a spherical domain with equally spaced nodes.

        The idea is simple. The requested number of equally spaced shells defines
        the linear spacing between nodes, dl. Use that linear spacing to fill the
        sphere, by filling up slices first, then stacks, then shells.
        
        The nodes are places in the _center_ of the shells, so that if rMin and
        rMax specify the boundaries of the domain, then a node at rMin+dl/2 is the
        center of the shell extending from rMin to rMin+dl. And similarly, a node
        placed at rMax-dl/2 is the center of the shell extending from rMax-dl to
        rMax.

        The stacks are built by placing nodes near the "north" and "south" poles,
        and as many equally spaced nodes as needed in between such that the
        linear distance between nodes on the great circle is close to dl. These
        colatitudes mark the _center_ of the stack.

        Each stack is then divided to equally spaced slices. But the slices don't
        all begin at a "prime meridian." Instead each is shifted a random angle.
        """

        loc_x = []
        loc_y = []
        loc_z = []
        loc_V = []
        loc_m = []
        loc_H = []

        self.linear_spacing = (self.rMax-self.rMin)/(self.nLayers) 
        h0 = 1.0/(self.linear_spacing*self.nNodePerh)
        nominalH = SymTensor3d(h0,  0.0, 0.0,
                               0.0, h0,  0.0,
                               0.0, 0.0, h0)

        # Fill up shells...
        dr = self.linear_spacing
        shells = np.linspace(self.rMin+dr/2, self.rMax-dr/2, self.nLayers)
        for r in shells:
            # With stacks...
            dG = dr/r
            nGs = max(int(pi/dG), 2)
            stacks = np.linspace(0.0+dG/2, pi-dG/2 ,nGs)
            for G in stacks:
                # Made of slices.
                dq = dr/(r*sin(G))
                nqs = max(int(2*pi/dq), 2)
                dq = 2*pi/nqs # A little bootstrapping for the smallest stacks...
                slices = np.linspace(0.0+dq/2, 2*pi-dq/2, nqs)
                slices += np.random.uniform(0,pi/4)
                for q in slices:
                    r0, r1 = r - dr/2, r + dr/2
                    G0, G1 = G - dG/2, G + dG/2
                    q0, q1 = q - dq/2, q + dq/2
                    dV = (1./3 * (r1**3-r0**3) * (cos(G0)-cos(G1)) * (q1-q0))
                    loc_x.append(r*sin(G)*cos(q))
                    loc_y.append(r*sin(G)*sin(q))
                    loc_z.append(r*cos(G))
                    loc_V.append(dV)
                    loc_m.append(self.rho * dV)
                    loc_H.append(nominalH)
                    pass
                pass
            pass

        # Save into the object's fields.
        self.x = loc_x
        self.y = loc_y
        self.z = loc_z
        self.V = loc_V
        self.m = loc_m
        self.H = loc_H
        
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
                 scale = 1.0,
                 rMin = 0.0,
                 rMax = 1e200,
                 nNodePerh = 2.01,
                 EOS = None,):
        """Class constructor for the HCP lattice node generator. 
        
        The generated coordinates are those of the CENTER of the node. Also note
        that because the lines and layers of the HCP lattice do not end at equal
        coordinate limits, the volume enclosed by the lattice is only 
        approximately equal to scale**3. 
           
        Parameters
        ----------
        nx : int > 0
            Number of nodes across domain. There are no corresponding ny or nz
            because we want to ensure a cubic lattice with equally spaced nodes
            in all directions.
        rho : float > 0
            Density used to assign node masses. For now, this is a constant.
        scale : float > 0, optional
            Linear scale of lattice. The lattice will be built with one corner
            placed at (0,0,0) and the opposite corner at (scale,scale,scale),
            and then translated to put the lattice center of mass at (0,0,0).
            This makes it simpler for a user to place the lattice in her 
            domain. Default is 1.0.
        rMin : float >=0, optional
            After lattice is built, nodes whose distance from the center of the
            lattice is less than rMin will be culled. Default is 0.0.
        rMax : float > 0, optional
            After lattice is built, nodes whose distance from the center of the
            lattice is greater than or equal to rMax will be culled. Default is 
            1e200.
        nNodePerh : float > 1.0, optional
            Nodes are assigned an inverse smoothing length of 1/(d*nNodePerh),
            where d is the lattice spacing. Default is 2.01.
        EOS : Spheral.SolidMaterial.SolidEquationOfState, optional.
            Place holder, for future use. Default is None.
        """

        # Some assertions for convenience. Not supposed to be an airtight seal.
        assert isinstance(nx,int) and nx > 1
        assert isinstance(rho,float) and rho > 0.0
        assert isinstance(scale,float) and scale > 0.0
        assert isinstance(rMin,float) and isinstance(rMax,float)
        assert rMax > rMin >= 0.0
        assert isinstance(nNodePerh,float) and nNodePerh >= 1.0

        # Store key parameters in the generator object.
        self.nx = nx
        self.rho = rho # will become a list when we get the node count.
        self.scale = scale
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

        # Convenience and diagnostic fields. (Don't overdo.)
        self.lattice_spacing = scale/(self.nx)
        self.lattice_volume = scale**3
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
        """Generate hcp lattice positions. 

        Each node will have 12 nearest neighbors, all a distance d from it, where
        d, the lattice spacing, is the diameter of the spheres that can be placed
        on the lattice in closest packing.
        """

        loc_x = []
        loc_y = []
        loc_z = []
        loc_V = []
        loc_m = []
        loc_H = []

        d = self.lattice_spacing
        r = d/2
        xstep = d
        ystep = sqrt(3)/2 * d
        zstep = sqrt(6)/3 * d
        nx = self.nx
        ny = int(nx * xstep/ystep) + 1
        nz = int(nx * xstep/zstep) + 1
        
        nominalCellVolume = self.lattice_volume/(nx*ny*nz)
        nominalCellMass = nominalCellVolume * self.rho
        h0 = 1.0/(d*self.nNodePerh)
        nominalH = SymTensor3d(h0,  0.0, 0.0,
                               0.0, h0,  0.0,
                               0.0, 0.0, h0)
        
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

                    loc_x.append(x)
                    loc_y.append(y)
                    loc_z.append(z)
                    loc_V.append(nominalCellVolume)
                    loc_m.append(nominalCellMass)
                    loc_H.append(nominalH)
                    pass
                pass
            pass

        # Translate the lattice to put the CoM at the origin. We do this in a
        # separate step to keep the original HCP calculation cleaner.
        xCM = sum(loc_x)/len(loc_x)
        yCM = sum(loc_y)/len(loc_y)
        zCM = sum(loc_z)/len(loc_z)
        loc_x = [x-xCM for x in loc_x]
        loc_y = [y-yCM for y in loc_y]
        loc_z = [z-zCM for z in loc_z]

        # Finally, chisel away a spherical shell.
        for k in range(len(loc_x)):
            R = hypot(loc_x[k], hypot(loc_y[k], loc_z[k]))
            if self.rMin <= R < self.rMax:
                self.x.append(loc_x[k])
                self.y.append(loc_y[k])
                self.z.append(loc_z[k])
                self.V.append(loc_V[k])
                self.m.append(loc_m[k])
                self.H.append(loc_H[k])
                pass
            pass
        
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
