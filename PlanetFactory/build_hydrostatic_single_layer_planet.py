#! /proj/nmovshov_hindmost/collisions/SPHERAL/bin/python
#-------------------------------------------------------------------------------
# Set up a single-layer, fluid planet and run to hydrostatic equilibrium.
#
# This script serves as a template for equilibrating a fluid, single material
# planet, with a specified mass and radius.
# Copy this file to a separate folder and modify with choices of planet mass,
# radius, and material. The parameter nxPlanet controls the run resolution. You
# will need to consider the expected time scale to run to, and you may need to
# tweak the cooldown frequency and strength as the planet approaches equilibrium.
#
# To run as an executable script, check that the shebang line points to the full
# path to spheral's python.
#-------------------------------------------------------------------------------
from math import *
import sys, os, shutil
import random
import mpi # Mike's simplified mpi wrapper
from SolidSpheral3d import *
from findLastRestart import findLastRestart
from VoronoiDistributeNodes import distributeNodes3d
from NodeHistory import NodeHistory
from GenerateNodeDistribution3d import GenerateNodeDistribution3d
pcsbase = '' # Edit this with full path to <pcs> if you see an ImportError.
sys.path += ['..',pcsbase,os.getenv('PCSBASE','')]
import shelpers # My module of some helper functions
import PlanetNodeGenerators # New experimental node generators

#-------------------------------------------------------------------------------
# NAV SETUP
# Physical parameters are in MKS units unless otherwise specified.
#-------------------------------------------------------------------------------

# Job name and description
jobName = 'onelayer'
jobDesc = "Hydrostatic equilibrium of a single-layer, fluid planet."
print '\n', jobName.upper(), '-', jobDesc.upper()

# Planet parameters
rPlanet = 1000e3             # Initial guess for planet radius (m)
rhoPlanet = 920.             # Initial guess for planet density (kg/m^3)
matPlanet = 'h2oice'         # Planet material (see <pcs>/MATERIALS.md for options)
mPlanet = rhoPlanet*4.0*pi/3.0*rPlanet**3
gravTime = 1/sqrt(MKS().G*rhoPlanet)

# Cooldown mechanism
cooldownMethod = 'dashpot'   # 'dashpot' or 'stomp' 
cooldownPower = 0.1          # Dimensionless cooldown "strength" >=0
cooldownFrequency = 1        # Cycles between application (use 1 with dashpot)
                             # * With 'stomp' method, 0<=power<=1

# Times, simulation control, and output
nxPlanet = 40                # Nodes across diameter of planet (run "resolution")
steps = None                 # None or number of steps to advance (overrides time)
goalTime = 1*gravTime        # Time to advance to (sec)
dtInit = 0.02                # Initial guess for time step (sec)
vizTime = 0.2*goalTime       # Time frequency for dropping viz files (sec)
vizCycle = None              # Cycle frequency for dropping viz files
outTime = vizTime            # Time between running output routine (sec)
outCycle = None              # Cycles between running output routine

# Node list parameters
nPerh = 2.01                 # Nominal number of nodes per smoothing scale
hmin = 1.0                   # Minimum smoothing length (fraction of nominal)
hmax = 1.0                   # Maximum smoothing length (fraction of nominal)
rhomax = 1e+1*rhoPlanet      # Upper bound on node density (kg/m^3)
generator_type = 'hcp'       # Node generator to use. 'hcp'|'old'|'shells'
hmin *= nPerh*2*rPlanet/nxPlanet
hmax *= nPerh*2*rPlanet/nxPlanet
rhomin = mPlanet/nxPlanet**3/hmax**3
universeEdge = 2*rPlanet

# Gravity parameters
softLength = 1.0             # Gravity softening length (fraction of nominal H)
opening = 1.0                # Opening parameter for gravity tree walk
fdt = 0.1                    # Time step multiplier (dt=fdt*sqrt(softlength/a))
G = MKS().G
softLength *= nPerh*2*rPlanet/nxPlanet

# More simulation parameters
dtGrowth = 2.0               # Maximum growth factor for time step per cycle 
dtMin = 0                    # Minimum allowed time step (sec)
dtMax = 0.1*goalTime         # Maximum allowed time step (sec)
verbosedt = False            # Verbose reporting of the time step criteria 
maxSteps = 2400              # Maximum allowed steps for simulation advance
statsStep = None             # Frequency for sampling conservation statistics etc.
redistributeStep = 8000      # Frequency to load balance problem from scratch
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
if cooldownFrequency is not None:
    assert 0 <= cooldownPower, "cool DOWN not up"
    if cooldownMethod is 'stomp':
        assert 0 <= cooldownPower <= 1.0, "stomp fraction is 0-1"
    assert type(cooldownFrequency) is int and cooldownFrequency > 0, "very funny"
    assert cooldownMethod in ['dashpot','stomp'], "unknown cooldown method"
    assert (cooldownFrequency == 1) or (not(cooldownMethod is 'dashpot')),\
            "dashpot cooling method requires frequency=1"
assert (outTime is None) or (outCycle is None),\
        "output on both time and cycle is confusing"
assert generator_type in ['hcp', 'shells', 'old']

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
hminratio = 1.0
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
# eos is determined by the  material string (see <pcs>MATERIALS.md for options).
#-------------------------------------------------------------------------------
eosPlanet = None

# Most eos constructors need to know about units. We usually use MKS.
units = PhysicalConstants(1.0, # unit length in meters
                          1.0, # unit mass in kilograms
                          1.0) # unit time in seconds

# Construct and verify planet eos
eosPlanet = shelpers.construct_eos_for_material(matPlanet,units)
assert eosPlanet is not None
assert eosPlanet.valid()

# Optionally, provide non-default values to the following
#eosPlanet.etamin_solid = 0.94 # default is 0.94
#eosPlanet.minimumPressure = 0.0 # default is -1e+200

#-------------------------------------------------------------------------------
# NAV Restarts and output directories
# Here we create the output directories, and deal with restarted runs if any.
#-------------------------------------------------------------------------------
# Name directories and files.
jobDir = os.path.join(baseDir, 
                       'mPlanet=%0.2g' % mPlanet,
                       'eosPlanet=%d' % eosPlanet.uid,
                       'nxPlanet=%i' % nxPlanet,
                       'np=%i' % mpi.procs,
                       )
restartDir = os.path.join(jobDir, 'restarts', 'proc-%04i' % mpi.rank)
vizDir = os.path.join(jobDir, 'viz')
outDir = os.path.join(jobDir, 'output')
logDir = os.path.join(jobDir, 'logs')
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
    if not os.path.exists(logDir):
        os.makedirs(logDir)
mpi.barrier()
if not os.path.exists(restartDir):
    os.makedirs(restartDir)
mpi.barrier()

# If we're restarting, find the set of most recent restart files.
if restoreCycle is None:
    restoreCycle = findLastRestart(restartName)

# Here's a quick way to save a record of parameters used in this run.
shutil.copyfile(__file__,logDir+'/{}.ini.{}'.format(jobName,restoreCycle))

#-------------------------------------------------------------------------------
# NAV Node construction
# Here we create and populate a node list with initial conditions. In spheral, the
# construction order is as follows:
# 1. Create an empty node list with fields matching the size and type of problem.
# 2. Create a "generator" that will decide what values to give all field variables
#    of node i. Normally we start with one of the simple, stock generators, and
#    modify the x,y,z,vx,vy,vz,rho,U values to suit our initial conditions.
# 3. Distribute, using the (nodeList, generator) pair, among ranks. The generator
#    will be used to fill values in the node list, and then discarded. 
#-------------------------------------------------------------------------------
# Create the node list.
planet = makeFluidNodeList('planet', eosPlanet, 
                           nPerh = nPerh, 
                           xmin = -universeEdge*Vector.one, # (probably unnecessary)
                           xmax =  universeEdge*Vector.one, # (probably unnecessary)
                           hmin = hmin,
                           hmax = hmax,
                           rhoMin = rhomin,
                           rhoMax = rhomax,
                           hminratio = hminratio,
                           topGridCellSize = universeEdge
                           )
planet.eos_id = eosPlanet.uid
nodeSet = [planet]

# Unless restarting, create the generator and set initial field values.
if restoreCycle is None:
    # Create a basic, usually constant density generator.
    if generator_type == 'old':
        planetGenerator = GenerateNodeDistribution3d(nxPlanet, nxPlanet, nxPlanet,
                            rhoPlanet,
                            distributionType = 'lattice',
                            xmin = (-rPlanet, -rPlanet, -rPlanet),
                            xmax = ( rPlanet,  rPlanet,  rPlanet),
                            rmin = 0.0,
                            rmax = rPlanet,
                            nNodePerh = nPerh)
        for k in range(planetGenerator.localNumNodes()):
            planetGenerator.x[k] *= 1.0 + random.uniform(-0.02, 0.02)
            planetGenerator.y[k] *= 1.0 + random.uniform(-0.02, 0.02)
            planetGenerator.z[k] *= 1.0 + random.uniform(-0.02, 0.02)
            pass
        pass
    elif generator_type == 'hcp':
        planetGenerator = PlanetNodeGenerators.HexagonalClosePacking(
                            nx = nxPlanet,
                            rho = rhoPlanet,
                            scale = 2*rPlanet,
                            rMin = 0.0,
                            rMax = rPlanet,
                            nNodePerh = nPerh)
        pass
    elif generator_type == 'shells':
        planetGenerator = PlanetNodeGenerators.EqualSpacingSphericalShells(
                            nLayers = nxPlanet/2,
                            rho = rhoPlanet,
                            rMin = 0.0,
                            rMax = rPlanet,
                            nNodePerh = nPerh)
        pass
    else:
        print "unknown generator type"
        sys.exit(1)
        pass

    # Tweak density profile to start closer to equilibrium.
    planetGenerator.EOS = eosPlanet
    shelpers.hydrostaticize_one_layer_planet(planetGenerator)
    
    # Fill node list using generator and distribute to ranks.
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
    
    pass # end restoreCycle branching

# The spheral controller needs a DataBase object to hold the node lists.
db = DataBase()
for n in nodeSet:
    db.appendNodeList(n)
del n

#-------------------------------------------------------------------------------
# NAV Spheral's simulation structure
# Here we construct the objects that make up spheral's simulation hierarchy.
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

# Create the kernel function for SPH.
WT = TableKernel(BSplineKernel(), 1000)

# Create the artificial viscosity object.
q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
q.epsilon2 = epsilon2
q.negligibleSoundSpeed = negligibleSoundSpeed
q.csMultiplier = csMultiplier

# Create the hydro package.
hydro = HydroConstructor(W = WT,
                         Q = q,
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
integrator = CheapSynchronousRK2Integrator(db)
integrator.appendPhysicsPackage(gravity)
integrator.appendPhysicsPackage(hydro)
integrator.lastDt = dtInit
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
#  * cooldown() - slow and cool internal nodes [cycle based]
#  * mOutput() - a generic access routine, usually a pickle of node list or some
#               calculated value of interest [cycle or time based]
#-------------------------------------------------------------------------------
def mOutput(stepsSoFar,timeNow,dt):
    mFileName="{0}-{1:05d}-{2:g}.{3}".format(
              jobName, stepsSoFar, timeNow, 'fnl.gz')
    shelpers.pflatten_node_list_list(nodeSet, outDir + '/' + mFileName)
    pass
if not outCycle is None:
    control.appendPeriodicWork(mOutput,outCycle)
if not outTime is None:
    control.appendPeriodicTimeWork(mOutput,outTime)

def cooldown(stepsSoFar,timeNow,dt):
    nbGlobalNodes = mpi.allreduce(sum([nl.numInternalNodes for nl in nodeSet]),
                                  mpi.SUM)
    massScale = mPlanet/nbGlobalNodes
    timeScale = 0.1*gravTime
    dashpotParameter = cooldownPower*massScale/timeScale
    for nl in nodeSet:
        v = nl.velocity()
        m = nl.mass()
        u = nl.specificThermalEnergy()
        if cooldownMethod == 'dashpot':
            for k in range(nl.numInternalNodes):
                v[k] *= 1.0 - min(dashpotParameter*dt/m[k], 1)
                u[k] *= 0.0 #TODO: maybe improve this
                pass
            pass
        elif cooldownMethod == 'stomp':
            for k in range(nl.numInternalNodes):
                v[k] *= 1.0 - cooldownPower
                u[k] *= 0.0 #TODO maybe improve this
                pass
            pass
        pass
    pass
control.appendPeriodicWork(cooldown,cooldownFrequency)

#-------------------------------------------------------------------------------
# NAV Launch simulation
# The simulation can be run for a specified number of steps, or a specified time
# in seconds.
#-------------------------------------------------------------------------------
# Save initial state in a flattened node list (.fnl) file.
mOutput(control.totalSteps, control.time(), control.lastDt())

# Report domain distribution statistics on first run.
if restoreCycle is None:
    dbase = control.integrator.dataBase()
    print "Running with {} internal nodes and {} ghost nodes".format(
            mpi.allreduce(dbase.numInternalNodes, mpi.SUM),
            mpi.allreduce(dbase.numGhostNodes, mpi.SUM))

# And go.
if not steps is None:
    print "Advancing {} steps.".format(steps)
    control.step(steps)
    control.dropRestartFile()
    control.dropViz()
else:
    print "Running to {} seconds.".format(goalTime)
    control.advance(goalTime, maxSteps)
    control.dropRestartFile()
    control.dropViz()

#-------------------------------------------------------------------------------
# NAV Post processing tasks
# Here we can include tasks that will happen once, if and when the run is completed
# successfully. Things like saving flattened node lists and/or computed quantities.
#-------------------------------------------------------------------------------
# Save final state in a flattened node list (.fnl) file.
mOutput(control.totalSteps, control.time(), control.lastDt())

# Print current planet's (approximate) vitals.
#TODO

#-------------------------------------------------------------------------------
# NAV Final thoughts
# Here we may print a message if desired, or do any final action.
#-------------------------------------------------------------------------------
print "\n", jobName.upper(), "completed."
