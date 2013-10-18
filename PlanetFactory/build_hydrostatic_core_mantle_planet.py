#! /proj/nmovshov_hindmost/collisions/SPHERAL/bin/python
#-------------------------------------------------------------------------------
# Set up a two-layer, fluid planet and run to hydrostatic equilibrium.
#
# This script serves as a template for equilibrating a fluid, two-layer  planet,
# specifying core and mantle radius and density.
# Copy this file to a separate folder and modify with choices of layer densities,
# radii, and materials. The parameter nxPlanet controls the run resolution. You
# will need to consider the expected time scale to run to, and you may need to
# tweak the cooldown frequency and strength as the planet approaches equilibrium.
#
# To run as an executable script, check that the shebang line points to the full
# path to spheral's python.
#-------------------------------------------------------------------------------
from math import *
import sys, os
import random
import mpi # Mike's simplified mpi wrapper
import shelpers # My module of some helper functions
from SolidSpheral3d import *
from findLastRestart import findLastRestart
from VoronoiDistributeNodes import distributeNodes3d
from NodeHistory import NodeHistory
from GenerateNodeDistribution3d import GenerateNodeDistribution3d

#-------------------------------------------------------------------------------
# NAV SETUP
# Physical parameters are in MKS units unless otherwise specified.
#-------------------------------------------------------------------------------

# Job name and description
jobName = 'coremantle'
jobDesc = "Hydrostatic equilibrium of a two-layer, fluid planet."
print '\n', jobName.upper(), '-', jobDesc.upper()

# Planet parameters
rPlanet = 1000e3             # Initial guess for outer planet radius (m)
rCore = 500e3                # Initial guess for core radius (m)
matMantle = 'h2oice'         # Mantle material (see <uss>/MATERIALS.md for options)
rhoMantle = 900.             # Initial guess for mantle density (kg/m^3)
matCore = 'granite'          # Core material (see <uss>/MATERIALS.md for options)
rhoCore = 2680.              # Initial guess for core density (kg/m^3)
mPlanet = (4.0*pi/3.0) * (rhoCore*rCore**3 + rhoMantle*(rPlanet**3-rCore**3))
rhoPlanet = 3.0*mPlanet/(4.0*pi*rPlanet**3)
gravTime = 1/sqrt(MKS().G*rhoPlanet)

# Cooldown mechanism
cooldownMethod = 'dashpot'   # 'dashpot' or 'stomp' 
cooldownPower = 1.0          # Dimensionless cooldown "strength" >=0
cooldownFrequency = 1        # Cycles between application (use 1 with dashpot)
                             # * With 'stomp' method, 0<=power<=1

# Times, simulation control, and output
nxPlanet = 20                # Nodes across diameter of planet (run "resolution")
steps = None                 # None or advance a number of steps rather than to a time
goalTime = 20*gravTime       # Time to advance to (sec)
dtInit = 0.02                # Initial guess for time step (sec)
vizTime = 0.4*gravTime       # Time frequency for dropping viz files (sec)
vizCycle = None              # Cycle frequency for dropping viz files
outTime = vizTime            # Time between running output routine (sec)
outCycle = None              # Cycles between running output routine

# Node list parameters
nPerh = 1.51                 # Nominal number of nodes per smoothing scale
hmin = 1e-6*rPlanet          # Lower bound on smoothing length
hmax = 2e+0*rPlanet          # Upper bound on smoothing length
rhomin = 1e-6*rhoPlanet      # Lower bound on node density
rhomax = 1e+1*rhoPlanet      # Upper bound on node density

# Gravity parameters
softLength = 1.0/nxPlanet    # Fraction of planet radius as softening length
opening = 1.0                # Dimensionless opening parameter for gravity tree walk
fdt = 0.1                    # Time step multiplier (dt=fdt*sqrt(softlength/a))
softLength *= rPlanet
G = MKS().G

# More simulation parameters
dtGrowth = 2.0               # Maximum growth factor for time step per cycle 
dtMin = 0                    # Minimum allowed time step (sec)
dtMax = 0.1*goalTime         # Maximum allowed time step (sec)
verbosedt = False            # Verbose reporting of the time step criteria per cycle
maxSteps = 800               # Maximum allowed steps for simulation advance
statsStep = None             # Frequency for sampling conservation statistics and such
redistributeStep = 2000      # Frequency to load balance problem from scratch
restartStep = 200            # Frequency to drop restart files
restoreCycle = None          # If None, latest available restart cycle is selected
baseDir = jobName            # Base name for directory to store output in

#-------------------------------------------------------------------------------
# NAV Assertions
# This is a good place for a quick abort if some bad parameter choices are going
# to cause trouble later, in confusing ways. We assume that spheral constructors
# use their own assertions, so here we can validate just our own stuff. Another
# valid option would be to simply not worry about it, and let exceptions happen.
#-------------------------------------------------------------------------------
assert 0 <= cooldownPower, "cool DOWN not up"
if cooldownMethod is 'stomp':
    assert 0 <= cooldownPower <= 1.0, "stomp fraction is 0-1"
assert type(cooldownFrequency) is int and cooldownFrequency > 0, "very funny"
assert cooldownMethod in ['dashpot','stomp'], "unknown cooldown method"
assert (cooldownFrequency == 1) or (not(cooldownMethod is 'dashpot')),\
        "dashpot cooling method requires frequency=1"
assert (outTime is None) or (outCycle is None),\
        "output on both time and cycle is confusing"
assert rPlanet > rCore, "core means it's inside"
assert matMantle != matCore, "why not use the single material script then?"

#-------------------------------------------------------------------------------
# NAV Spheral hydro solver options
# These options for spheral's hydro mechanism are normally left alone.
#-------------------------------------------------------------------------------
HydroConstructor = ASPHHydro
Qconstructor = MonaghanGingoldViscosity
Cl = 1.0
Cq = 1.0
Qlimiter = False
balsaraCorrection = False
epsilon2 = 1e-2
negligibleSoundSpeed = 1e-4 # kind of arbitrary.
csMultiplier = 1e-4
hminratio = 0.1
limitIdealH = False
cfl = 0.5
useVelocityMagnitudeForDt = False
XSPH = True
epsilonTensile = 0.0
nTensile = 4
HEvolution = IdealH
densityUpdate = IntegrateDensity # (Sum|Integrate)Density
compatibleEnergyEvolution = True
rigorousBoundaries = False

#-------------------------------------------------------------------------------
# NAV Equation of state
# Here we construct equation-of-state objects, one per node list. The choice of 
# eos is determined by the  material string (see <uss>MATERIALS.md for options).
#-------------------------------------------------------------------------------
eosCore, eosMantle = None, None

# Most eos constructors need to know about units. We usually use MKS.
units = PhysicalConstants(1.0, # unit length in meters
                          1.0, # unit mass in kilograms
                          1.0) # unit time in seconds

# Construct and verify core eos
eosCore = shelpers.construct_eos_for_material(matCore,units)
assert eosCore is not None
assert eosCore.valid()

# Construct and verify mantle eos
eosMantle = shelpers.construct_eos_for_material(matMantle,units)
assert eosMantle is not None
assert eosMantle.valid()

# Optionally, provide non-default values to the following
eosCore.etamin = 0.94
eosCore.minimumPressure = -1e200
eosMantle.etamin = 0.94
eosMantle.minimumPressure = -1e200

#-------------------------------------------------------------------------------
# NAV Restarts and output directories
# Here we create the output directories, and deal with restarted runs if any.
#-------------------------------------------------------------------------------
# Name directories and files.
jobDir = os.path.join(baseDir, 
                       'rPlanet=%0.2g' % rPlanet,
                       'rCore=%0.2g' % rCore,
                       'eosCore=%d' % eosCore.uid,
                       'eosMantle=%d' % eosMantle.uid,
                       'nxPlanet=%i' % nxPlanet,
                       'np=%i' % mpi.procs,
                       )
restartDir = os.path.join(jobDir, 'restarts', 'proc-%04i' % mpi.rank)
vizDir = os.path.join(jobDir, 'viz')
outDir = os.path.join(jobDir, 'output')
restartName = os.path.join(restartDir, jobName)

# Check if the necessary directories exist.  If not, create them.
if mpi.rank == 0:
    if not os.path.exists(jobDir):
        os.makedirs(jobDir)
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
    if not os.path.exists(outDir):
        os.makedirs(outDir)
mpi.barrier()
if not os.path.exists(restartDir):
    os.makedirs(restartDir)
mpi.barrier()

# If we're restarting, find the set of most recent restart files.
if restoreCycle is None:
    restoreCycle = findLastRestart(restartName)

#-------------------------------------------------------------------------------
# NAV Node construction
# Here we create and populate node lists with initial conditions. In spheral, the
# construction order is as follows:
# 1. Create an empty node list with fields matching the size and type of problem.
# 2. Create a "generator" that will decide what values to give all field variables
#    of node i. Normally we start with one of the simple, stock generators, and
#    modify the x,y,z,vx,vy,vz,rho,U values to suit our initial conditions.
# 3. Distribute, using the (nodeList, generator) pair, among ranks. The generator
#    will be used to fill values in the node list, and then discarded. 
#-------------------------------------------------------------------------------
# Create the node lists.
core   = makeFluidNodeList('core', eosCore, 
                           nPerh = nPerh, 
                           xmin = -10.0*rCore*Vector.one, # (probably unnecessary)
                           xmax =  10.0*rCore*Vector.one, # (probably unnecessary)
                           hmin = hmin,
                           hmax = hmax,
                           rhoMin = rhomin,
                           rhoMax = rhomax,
                           )
core.eos_id = eosCore.uid

mantle = makeFluidNodeList('mantle', eosMantle, 
                           nPerh = nPerh, 
                           xmin = -10.0*rPlanet*Vector.one, # (probably unnecessary)
                           xmax =  10.0*rPlanet*Vector.one, # (probably unnecessary)
                           hmin = hmin,
                           hmax = hmax,
                           rhoMin = rhomin,
                           rhoMax = rhomax,
                           )
mantle.eos_id = eosMantle.uid

nodeSet = [core, mantle]

# Unless restarting, create the generator and set initial field values.
if restoreCycle is None:
    # Start with the stock generator.
    nxCore = int(nxPlanet*(rCore/rPlanet))
    coreGenerator   = GenerateNodeDistribution3d(nxCore, nxCore, nxCore,
                                          rhoCore,
                                          distributionType = 'lattice',
                                          xmin = (-rCore, -rCore, -rCore),
                                          xmax = ( rCore,  rCore,  rCore),
                                          rmin = 0.0,
                                          rmax = rCore,
                                          nNodePerh = nPerh)
    mantleGenerator = GenerateNodeDistribution3d(nxPlanet, nxPlanet, nxPlanet,
                                          rhoMantle,
                                          distributionType = 'lattice',
                                          xmin = (-rPlanet, -rPlanet, -rPlanet),
                                          xmax = ( rPlanet,  rPlanet,  rPlanet),
                                          rmin = rCore,
                                          rmax = rPlanet,
                                          nNodePerh = nPerh)

    # We disturb the lattice symmetry to avoid artificial singularities.
    for k in range(coreGenerator.localNumNodes()):
        coreGenerator.x[k] *= 1.0 + random.uniform(-0.02, 0.02)
        coreGenerator.y[k] *= 1.0 + random.uniform(-0.02, 0.02)
        coreGenerator.z[k] *= 1.0 + random.uniform(-0.02, 0.02)
        pass
    for k in range(mantleGenerator.localNumNodes()):
        mantleGenerator.x[k] *= 1.0 + random.uniform(-0.02, 0.02)
        mantleGenerator.y[k] *= 1.0 + random.uniform(-0.02, 0.02)
        mantleGenerator.z[k] *= 1.0 + random.uniform(-0.02, 0.02)
        pass

    # Fill node lists using generators and distribute to ranks.
    print "Starting node distribution..."
    distributeNodes3d((core, coreGenerator),
                      (mantle, mantleGenerator))
    nGlobalNodes = 0
    for n in nodeSet:
        print "Generator info for %s" % n.name
        print "   Minimum number of nodes per domain : ", \
              mpi.allreduce(n.numInternalNodes, mpi.MIN)
        print "   Maximum number of nodes per domain : ", \
              mpi.allreduce(n.numInternalNodes, mpi.MAX)
        print "               Global number of nodes : ", \
              mpi.allreduce(n.numInternalNodes, mpi.SUM)
        nGlobalNodes += mpi.allreduce(n.numInternalNodes, mpi.SUM)
    del n
    print "Total number of (internal) nodes in simulation: ", nGlobalNodes
    
    pass
# The spheral controller needs a DataBase object to hold the node lists.
db = DataBase()
for n in nodeSet:
    db.appendNodeList(n)
del n
