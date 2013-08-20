#! /proj/nmovshov_hindmost/collisions/SPHERAL/bin/python
#-------------------------------------------------------------------------------
# Run a polytropic fluid planet to hydrostatic equilibrium.
#
# This script demonstrates the interaction of gravity and hydrodynamics inside a
# planet sized body. Using an index 1 polytropic equation-of-state, the planet
# should converge to a known density profile. But some tricks need to be employed
# to help the inherently dynamic spheral reach a more-or-less static equilibrium.
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
from findLastRestart import *
from VoronoiDistributeNodes import distributeNodes3d
from NodeHistory import NodeHistory
from GenerateNodeDistribution3d import GenerateNodeDistribution3d

#-------------------------------------------------------------------------------
# NAV SETUP
# Physical parameters are in MKS units unless otherwise specified.
#-------------------------------------------------------------------------------

# Job name and description
jobName = 'polystatic'
jobDesc = "Hydrostatic equilibrium of a polytropic planet."
print '\n', jobName.upper(), '-', jobDesc.upper()

# Planet properties
rPlanet = 11.2            # Initial guess for radius of planet (earth radii)
mPlanet = 318             # Mass of planet (earth masses)
polytrope_K  = 2e5        # Polytropic constant (varies)
polytrope_n  = 1          # Polytropic index (n=1)
polytrope_mu = 2.2e-3     # Mean molecular weight (kg/mole)
mPlanet *= 5.972e24
rPlanet *= 6371.0e3
rhoPlanet = 3.0*mPlanet/(4.0*pi*rPlanet**3)

# Times, simulation control, and output
steps = None              # None or advance a number of steps rather than to a time
goalTime = 24000          # Time to advance to (sec)
dt = 20                   # Initial guess for time step (sec)
vizTime = 1200            # Time frequency for dropping viz files (sec)
vizCycle = None           # Cycle frequency for dropping viz files

# Node seeding parameters ("resolution")
nxPlanet = 20             # Number of nodes across the diameter of the target
nPerh = 1.51              # Nominal number of nodes per smoothing scale
hmin = 1.0e-6*rPlanet     # Lower bound on smoothing length
hmax = 1.0e-1*rPlanet     # Upper bound on smoothing length
rhomin = 0.001*rhoPlanet  # Lower bound on node density
rhomax = 4.0*rhoPlanet    # Upper bound on node density

# Gravity parameters
softLength = 1.0e-5       # Fraction of planet radius as softening length
opening = 1.0             # Dimensionless opening parameter for gravity tree walk
fdt = 0.1                 # Gravity time step multiplier
softLength *= rPlanet
G = MKS().G

# More simulation parameters
dtGrowth = 2.0            # Maximum growth factor for time step in a cycle (dimensionless)
dtMin = 2                 # Minimum allowed time step (sec)
dtMax = 1000.0*dt         # Maximum allowed time step (sec)
verbosedt = False         # Verbose reporting of the time step criteria per cycle
maxSteps = 1000           # Maximum allowed steps for simulation advance
statsStep = None          # Frequency for sampling conservation statistics and such
redistributeStep = 2000   # Frequency to load balance problem from scratch
restartStep = 200         # Frequency to drop restart files
restoreCycle = None       # If None, latest available restart cycle is selected
baseDir = jobName         # Base name for directory to store output in

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
densityUpdate = RigorousSumDensity # Sum is best for fluids, integrate for solids.
compatibleEnergyEvolution = True
rigorousBoundaries = False

#-------------------------------------------------------------------------------
# NAV Equation of state
#-------------------------------------------------------------------------------
eosPlanet = PolytropicEquationOfStateMKS3d(polytrope_K,
                                           polytrope_n,
                                           polytrope_mu)
assert eosPlanet.valid(), "equation of state construction failed"

#-------------------------------------------------------------------------------
# NAV Restarts and output directories
# Here we create the output directories, and deal with restarted runs if any.
#-------------------------------------------------------------------------------
# Name directories and files.
jobDir = os.path.join(baseDir, 
                       'index=%g' % polytrope_n,
                       'const=%g' % polytrope_K,
                       'nxPlanet=%i' % nxPlanet,
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
# NAV Here we construct a node list for a spherical stationary planet
#-------------------------------------------------------------------------------
# Create the NodeList.
planet = makeFluidNodeList("planet", eosPlanet, 
                           nPerh = nPerh, 
                           xmin = -10.0*rPlanet*Vector.one,
                           xmax =  10.0*rPlanet*Vector.one,
                           hmin = hmin,
                           hmax = hmax,
                           rhoMin = rhomin,
                           rhoMax = rhomax,
                           )
nodeSet = [planet]

# Generate nodes
if restoreCycle is None:
    print "Generating node distribution."
    from GenerateNodeDistribution3d import GenerateNodeDistribution3d

    planetGenerator = GenerateNodeDistribution3d(nxPlanet, nxPlanet, nxPlanet,
                                                 rhoPlanet,
                                                 distributionType = "lattice",
                                                 xmin = (-rPlanet, -rPlanet, -rPlanet),
                                                 xmax = ( rPlanet,  rPlanet,  rPlanet),
                                                 rmax = rPlanet,
                                                 nNodePerh = nPerh)

# Disturb the symmetry with some random noise to avoid artificial waves
    for k in range(planetGenerator.localNumNodes()):
        planetGenerator.x[k] *= 1.0 + 0.04*random.random()
        planetGenerator.y[k] *= 1.0 + 0.04*random.random()
        planetGenerator.z[k] *= 1.0 + 0.04*random.random()

# Distribute nodes across ranks
    print "Starting node distribution..."
    distributeNodes3d((planet, planetGenerator))
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
    
# Construct a DataBase to hold our node list.
db = DataBase()
for n in nodeSet:
    db.appendNodeList(n)
del n

#-------------------------------------------------------------------------------
# NAV Here we create the various physics objects
#-------------------------------------------------------------------------------

# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
WT = TableKernel(BSplineKernel(), 1000)
WTPi = WT

# Create the gravity object
gravity = OctTreeGravity(G = G, 
                         softeningLength = softLength, 
                         opening = opening, 
                         ftimestep = fdt)

# Construct the artificial viscosity.
q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
q.epsilon2 = epsilon2
q.negligibleSoundSpeed = negligibleSoundSpeed
q.csMultiplier = csMultiplier

# Construct the hydro physics object.
hydro = HydroConstructor(WT,
                         WTPi,
                         q,
                         cfl = cfl,
                         useVelocityMagnitudeForDt = useVelocityMagnitudeForDt,
                         compatibleEnergyEvolution = compatibleEnergyEvolution,
                         gradhCorrection = False,
                         densityUpdate = densityUpdate,
                         HUpdate = HEvolution,
                         XSPH = XSPH,
                         epsTensile = epsilonTensile,
                         nTensile = nTensile)

# Construct a time integrator.
integrator = SynchronousRK2Integrator(db)
integrator.appendPhysicsPackage(gravity)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.verbose = verbosedt
integrator.rigorousBoundaries = rigorousBoundaries

# Build the controller.
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            redistributeStep = redistributeStep,
                            restartBaseName = restartName,
                            restoreCycle = restoreCycle,
                            vizBaseName = baseName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime)

#-------------------------------------------------------------------------------
# NAV MIDPROCESS Here we register optional work to be done mid-run
#-------------------------------------------------------------------------------
def midprocess(stepsSoFar,timeNow,dt):
    # stop and cool all nodes
    vref = planet.velocity()
    uref = planet.specificThermalEnergy()
    for k in range(planet.numNodes):
        vref[k] *= cooldownFactor
        uref[k] *= cooldownFactor
    pass
frequency=cooldownFrequency
control.appendPeriodicWork(midprocess,frequency)

#-------------------------------------------------------------------------------
# NAV Here we launch the simulation
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
    #raise ValueError, ("Completed %i steps." % steps)

else:
    control.advance(goalTime, maxSteps)
    control.dropRestartFile()
    #control.step() # One more step to ensure we get the final viz dump.

#-------------------------------------------------------------------------------
# NAV Here we do any post processing
#-------------------------------------------------------------------------------
shelpers.spickle_node_list(planet, jobName+'.dat')
#from IPython import embed
#if mpi.rank == 0:
#    embed() # uncomment to start an interactive session when the run completes

