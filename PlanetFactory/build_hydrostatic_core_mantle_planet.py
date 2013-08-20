#! /proj/nmovshov_hindmost/collisions/SPHERAL/bin/python
#-------------------------------------------------------------------------------
# Set up a two-layer, fluid planet and run to hydrostatic equilibrium.
#
# This script can serve as a template for equilibrating a fluid, two layer planet,
# with various core/mantle materials and size ratios.
# Copy this file to a separate folder and modify with your initial conditions. You
# will need to consider the expected time scale to run, and you may need to tweak
# the cooldown frequency and strength as the planet approaches equilibrium.
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
# Modify these variables to suit your specific problem. Physical parameters are
# in MKS units unless otherwise specified.
#-------------------------------------------------------------------------------

# Job name and description
jobName = 'coremantle'
jobDesc = "Hydrostatic equilibrium of a two-layer, fluid planet."
print '\n', jobName.upper(), '-', jobDesc.upper()

# Planet properties
rPlanet = 500e3           # Initial guess for outer radius of planet (m)
rCore   = 250e3           # Initial guess for radius of inner core (m)
matMantle = 'pure ice'    # Mantle material, see ../MATERIALS.md for options
matCore   = 'basalt'      # Core material, see ../MATERIALS.md for options
rhoMantle = 917           # Initial mantle mean density (kg/m^3)
rhoCore = 2700            # Initial core mean density (kg/m^3)
mPlanet = (4.0*pi/3.0) * (rhoCore*rCore**3 + rhoMantle*(rPlanet**3 - rCore**3))

# Times, simulation control, and output
steps = None              # None or advance a number of steps rather than to a time
goalTime = 1800           # Time to advance to (sec)
dt = 2                    # Initial guess for time step (sec)
vizTime = 200             # Time frequency for dropping viz files (sec)
vizCycle = None           # Cycle frequency for dropping viz files
cooldownFrequency = 1     # None or cycles between "cooldowns" (v=0, U=0)
cooldownFactor = 0.8      # 0.0-1.0 multiplier of velocity and energy during cooldown

# Node seeding parameters ("resolution")
nxPlanet = 20             # Number of nodes across the diameter of the target
nPerh = 1.51              # Nominal number of nodes per smoothing scale
hmin = 1.0e-6*rPlanet     # Lower bound on smoothing length
hmax = 1.0e-1*rPlanet     # Upper bound on smoothing length
rhomin = 0.01*rhoMantle   # Lower bound on node density
rhomax = 4.0*rhoCore      # Upper bound on node density

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
restartStep = 100         # Frequency to drop restart files
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
# NAV Equations of state
# Here we construct an eos object for each node list. In this case, one fore core,
# one for mantle. The choice of eos is determined by the material string. See
# ../MATERIALS.md for the available options.
# TODO: fix ANEOS, currently only tillotson works.
#-------------------------------------------------------------------------------
eosCore, eosMantle = None, None
# Most eos constructors take a units object, we usually use MKS.
units = PhysicalConstants(1.0,   # Unit length in meters
                          1.0,   # Unit mass in kg
                          1.0)   # Unit time in seconds

# Use tillotson EOS for many geologic materials.
mats = ['granite', 'basalt', 'nylon', 'pure ice', '30% silicate ice', 'water']
etamin, etamax = 0.01, 100.0
if matMantle.lower() in mats:
    eosMantle = TillotsonEquationOfState(matMantle,etamin,etamax,units)
if matCore.lower() in mats:
    eosCore = TillotsonEquationOfState(matCore,etamin,etamax,units)

# Verify valid EOSs (currently only tillotson works).
if eosCore is None or eosMantle is None:
    raise ValueError("invalid material selection for core and/or mantle")
if not (eosCore.valid() and eosMantle.valid()):
    raise ValueError("core and/or mantle eos construction failed")

#-------------------------------------------------------------------------------
# NAV Restarts and output directories
# Here we create the output directories, and deal with restarted runs if any.
#-------------------------------------------------------------------------------
# Name directories and files.
jobDir = os.path.join(baseDir, 
                       'core=%s' % matCore,
                       'mantle=%s' % matMantle,
                       'nxPlanet=%d' % nxPlanet,
                       )
restartDir = os.path.join(jobDir, 'restarts', 'proc-%04d' % mpi.rank)
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
# Here we create and populate the node lists for the core and mantle. In spheral,
# the construction order is as follows:
# 1. Create an empty node list with fields that match the size and type of problem.
# 2. Create a "generator" that will decide what values to give to all field variables
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

mantle = makeFluidNodeList('mantle', eosMantle, 
                           nPerh = nPerh, 
                           xmin = -10.0*rPlanet*Vector.one, # probably unnecessary
                           xmax =  10.0*rPlanet*Vector.one, # probably unnecessary
                           hmin = hmin,
                           hmax = hmax,
                           rhoMin = rhomin,
                           rhoMax = rhomax,
                           )
nodeSet = [core, mantle]

# Unless restarting, create the generators and set initial field values.
if restoreCycle is None:
    # Start with the stock generator.
    coreGenerator   = GenerateNodeDistribution3d(nxPlanet, nxPlanet, nxPlanet,
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
                                                 rmin = rCore + rPlanet/nxPlanet,
                                                 rmax = rPlanet,
                                                 nNodePerh = nPerh)
     
    # Modify geometry. 
    # We disturb the lattice symmetry to avoid artificial singularities.
    for k in range(coreGenerator.localNumNodes()):
        coreGenerator.x[k]   *= 1.0 + random.uniform(-0.02, 0.02)
        coreGenerator.y[k]   *= 1.0 + random.uniform(-0.02, 0.02)
        coreGenerator.z[k]   *= 1.0 + random.uniform(-0.02, 0.02)
    for k in range(mantleGenerator.localNumNodes()):
        mantleGenerator.x[k] *= 1.0 + random.uniform(-0.02, 0.02)
        mantleGenerator.y[k] *= 1.0 + random.uniform(-0.02, 0.02)
        mantleGenerator.z[k] *= 1.0 + random.uniform(-0.02, 0.02)
    pass

    # Modify density.
    pass

    # Modify velocity.
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
    
# The spheral controller needs a DataBase object to hold the node lists.
db = DataBase()
for n in nodeSet:
    db.appendNodeList(n)
del n

#-------------------------------------------------------------------------------
# NAV Spheral's simulation structure
# Here we construct the objects that compose spheral's simulation hierarchy.
# These are:
#  * One or more physics packages (hydro, gravity, strength, damage)
#  * A time integrator of some flavor (usually a Runge-Kutta 2)
#  * The simulation controller
#-------------------------------------------------------------------------------

# Create the gravity package.
gravity = OctTreeGravity(G = G, 
                         softeningLength = softLength, 
                         opening = opening, 
                         ftimestep = fdt)

# Create the kernel functions for SPH.
WT = TableKernel(BSplineKernel(), 1000) # one for normal hydro
WTPi = WT                               # one for artificial viscosity

# Create the artificial viscosity object.
q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
q.epsilon2 = epsilon2
q.negligibleSoundSpeed = negligibleSoundSpeed
q.csMultiplier = csMultiplier

# Create the hydro package.
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

# Create the time integrator and attach the physics packages to it.
integrator = SynchronousRK2Integrator(db)
integrator.appendPhysicsPackage(gravity)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dt
integrator.dtMin = dtMin
integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.verbose = verbosedt
integrator.rigorousBoundaries = rigorousBoundaries

# Create the controller.
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            redistributeStep = redistributeStep,
                            restartBaseName = restartName,
                            restoreCycle = restoreCycle,
                            vizBaseName = jobName,
                            vizDir = vizDir,
                            vizStep = vizCycle,
                            vizTime = vizTime)

#-------------------------------------------------------------------------------
# NAV Periodic, mid-run actions
# Here we register optional work to be done mid-run. Mid-run processes can be time
# or cycle based. Here we use:
#  * cooldown() - slow and cool all nodes (including ghost!) [cycle based]
#-------------------------------------------------------------------------------
def cooldown(stepsSoFar,timeNow,dt):
    # stop and cool all(!) nodes
    v_core = core.velocity()
    u_core = core.specificThermalEnergy()
    for k in range(core.numNodes):
        v_core[k] *= cooldownFactor
        u_core[k] *= cooldownFactor
    v_mantle = mantle.velocity()
    u_mantle = mantle.specificThermalEnergy()
    for k in range(mantle.numNodes):
        v_mantle[k] *= cooldownFactor
        u_mantle[k] *= cooldownFactor
    pass
    # end cooldown()

frequency=cooldownFrequency
control.appendPeriodicWork(cooldown,frequency)

#-------------------------------------------------------------------------------
# NAV Launch simulation
# The simulation can be run for a specified number of steps, or a specified time
# in seconds.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
else:
    control.advance(goalTime, maxSteps)
    control.dropRestartFile()
    control.dropViz()

#-------------------------------------------------------------------------------
# NAV Post processing tasks
# Here we can include tasks that will happen once, if and when the run is completed
# successfully. Things like saving flattened node lists and/or computed quantities.
#-------------------------------------------------------------------------------
pass
