#! /proj/nmovshov_hindmost/collisions/SPHERAL/bin/python
#-------------------------------------------------------------------------------
# Set up a two-layer, fluid planet and run to hydrostatic equilibrium.
#
#-------------------------------------------------------------------------------
from math import *
import sys
import random
import mpi # Mike's simplified mpi wrapper
import shelpers # My module of some helper functions
from SolidSpheral3d import *
from findLastRestart import *
from VoronoiDistributeNodes import distributeNodes3d
from NodeHistory import NodeHistory

#-------------------------------------------------------------------------------
# NAV SETUP
# Modify these variables to suit your specific problem. Physical parameters are
# in MKS units unless otherwise specified.
#-------------------------------------------------------------------------------

# Job name and description
jobName = "coremantle"
jobDesc = "Hydrostatic equilibrium of a two-layer, fluid planet."
print jobName, '-', jobDesc

# Planet properties
rPlanet = 500e3           # (m) initial guess for outer radius of planet
rCore   = 250e3           # (m) initial guess for radius of inner core
matMantle = 'h2oice'      # mantle material, see ../MATERIALS.md for options
matCore   = 'basalt'      # core material, see ../MATERIALS.md for options
rhoMantle = 917           # (kg/m^3) initial mantle mean density
rhoCore = 2700            # (kg/m^3) initial core mean density
mPlanet = (4.0/3.0*pi) * (rhoCore*rCore**3 + rhoMantle*(rPlanet**3 - rCore**3))

# Times, simulation control, and output
steps = None              # None or advance a number of steps rather than to a time
goalTime = 24000          # Time to advance to (sec)
dt = 20                   # Initial guess for time step (sec)
vizTime = 600             # Time frequency for dropping viz files (sec)
vizCycle = None           # Cycle frequency for dropping viz files
cooldownFrequency = 1     # None or cycles between "cooldowns" (v=0, U=0)
cooldownFactor = 0.8      # 0.0-1.0 multiplier of velocity and energy during cooldown

# Node seeding parameters ("resolution")
nxPlanet = 40             # Number of nodes across the diameter of the target
nPerh = 1.51              # Nominal number of nodes per smoothing scale
hmin = 1.0e-6*rPlanet     # Lower bound on smoothing length
hmax = 1.0e-1*rPlanet     # Upper bound on smoothing length
rhomin = 0.01*rhoPlanet   # Lower bound on node density
rhomax = 4.0*rhoPlanet    # Upper bound on node density

# Gravity parameters
softLength = 1.0e-5       # (fraction of planet radius) softening length
opening = 1.0             # (dimensionless) opening parameter for gravity tree walk
fdt = 0.1                 # (dimensionless) gravity timestep multiplier
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
restartStep = 100         # Frequency to drop restart files
restoreCycle = None       # If None, latest available restart cycle is selected
baseDir = jobName         # Base name for directory to store output in

#-------------------------------------------------------------------------------
# NAV Options for spheral's hydro mechanism (normally left alone)
#-------------------------------------------------------------------------------
HydroConstructor = ASPHHydro
Qconstructor = MonaghanGingoldViscosity
Cl = 1.0
Cq = 1.0
Qlimiter = False
balsaraCorrection = False
epsilon2 = 1e-2
negligibleSoundSpeed = 1e-4 #TODO make physics based
csMultiplier = 1e-4
hminratio = 0.1
limitIdealH = False
cfl = 0.5
useVelocityMagnitudeForDt = False
XSPH = True
epsilonTensile = 0.3
nTensile = 4
HEvolution = IdealH
densityUpdate = RigorousSumDensity # Sum is best for fluids, integrate for solids
compatibleEnergyEvolution = True
rigorousBoundaries = False

#-------------------------------------------------------------------------------
# NAV Build polytropic EOS object
#-------------------------------------------------------------------------------
eosPlanet = PolytropicEquationOfStateMKS3d(polytrope_K,
                                           polytrope_n,
                                           polytrope_mu)

#-------------------------------------------------------------------------------
# NAV Here we determine if, and deal with, restarted runs
#-------------------------------------------------------------------------------
# Restart and output files.
dataDir = os.path.join(baseDir, 
                       "index=%g" % polytrope_n,
                       "const=%g" % polytrope_K,
                       "nxPlanet=%i" % nxPlanet,
                       )
restartDir = os.path.join(dataDir, "restarts", "proc-%04i" % mpi.rank)
vizDir = os.path.join(dataDir, "viz")
baseName = jobName
restartName = os.path.join(restartDir, baseName)

# Check if the necessary output directories exist.  If not, create them.
import os, sys
if mpi.rank == 0:
    if not os.path.exists(dataDir):
        os.makedirs(dataDir)
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
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


