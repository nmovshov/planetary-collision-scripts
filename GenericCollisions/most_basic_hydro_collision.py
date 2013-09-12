#! /proj/nmovshov_hindmost/collisions/SPHERAL/bin/python
#-------------------------------------------------------------------------------
# Run two fluid spheres into each other.
# 
# This script launches the most basic type of collision we can imagine. Two
# spherical fluid objects of arbitrary size collide with some specified velocity
# and impact parameter. In spheral terms, the target and impactor are of type
# FluidNodeList, and the only physics package attached to the integrator is the
# hydro package. Although the target and impactor may select a solid material 
# equation-of-state, they do not possess any elastic strength, and no damage 
# model is attached. No gravity is at play either.
#
# Although this type of collision is not really physically relevant, it may 
# still serve as a sanity check on collisions with more complex physics, and
# also help identify what physical process is most dominant in a given impact.
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
jobName = 'mostbasic'
jobDesc = "Pure hydro collision of fluid, single material spheres."
print '\n', jobName.upper(), '-', jobDesc.upper()

# Target parameters
rTarget = 220.0            # Target radius (m)
mTarget = 3.36e10          # Target mass (kg)
matTarget = 'pure ice'     # Target material (see uss/MATERIALS.md for options)
rhoTarget = 3.0*mTarget/(4.0*pi*rTarget**3)

# Impactor parameters
rImpactor = 100.0          # Impactor radius (m)
mImpactor = 4.28e9         # Impactor mass (kg)
matImpactor = 'pure ice'   # Impactor material (see uss/MATERIALS.md for options)
rhoImpactor = 3.0*mImpactor/(4.0*pi*rImpactor**3)

# Collision parameters
vImpact = 10               # Impact velocity (m/s)
angleImpact = 30           # Impact angle to normal (degrees)

# Times, simulation control, and output
nxTarget = 20              # Nodes across diameter of target (run "resolution")
steps = 0               # None or advance a number of steps rather than to a time
goalTime = 6               # Time to advance to (sec)
dtInit = 0.02              # Initial guess for time step (sec)
vizTime = None             # Time frequency for dropping viz files (sec)
vizCycle = 1            # Cycle frequency for dropping viz files
outTime = None             # Time between running output routine (sec)
outCycle = 1            # Cycles between running output routine

# Node list parameters
nPerh = 1.51               # Nominal number of nodes per smoothing scale
hmin = 1e-6*rTarget        # Lower bound on smoothing length
hmax = 1e+1*rTarget        # Upper bound on smoothing length
rhomin = 1e-4*rhoTarget    # Lower bound on node density
rhomax = 1e+8*rhoTarget    # Upper bound on node density

# More simulation parameters
dtGrowth = 2.0             # Maximum growth factor for time step per cycle (dimensionless)
dtMin = 0                  # Minimum allowed time step (sec)
dtMax = 0.1*goalTime       # Maximum allowed time step (sec)
verbosedt = True           # Verbose reporting of the time step criteria per cycle
maxSteps = 1000            # Maximum allowed steps for simulation advance
statsStep = None           # Frequency for sampling conservation statistics and such
redistributeStep = 2000    # Frequency to load balance problem from scratch
restartStep = 200          # Frequency to drop restart files
restoreCycle = None        # If None, latest available restart cycle is selected
baseDir = jobName          # Base name for directory to store output in

#-------------------------------------------------------------------------------
# NAV Assertions
# This is a good place for a quick abort if some bad parameter choices are going
# to cause trouble later, in confusing ways. We assume that spheral constructors
# use their own assertions, so here we can validate just our own stuff. Another
# valid option would be to simply not worry about it, and let exceptions happen.
#-------------------------------------------------------------------------------
assert 0 <= angleImpact < 90, "give impact angle in first quadrant (in degrees)"
assert (outTime is None) or (outCycle is None),\
        "output on both time and cycle is confusing"

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
# Here we construct equation-of-state objects, one per node list. In this case,
# one each for the target and impactor. The choice of eos is determined by the 
# material string (see <uss>MATERIALS.md for available options).
#-------------------------------------------------------------------------------
eosTarget, eosImpactor = None, None
# Most eos constructors need to know about units. We usually use MKS.
units = PhysicalConstants(1.0, # unit length in meters
                          1.0, # unit mass in kilograms
                          1.0) # unit time in seconds

# Select and construct target eos
if matTarget.lower() in shelpers.material_strings['tillotson']:
    etamin, etamax = rhomin/rhoTarget, rhomax/rhoTarget
    pext, pmin, pmax = 0.0, 0.0, 1e200
    eosTarget = TillotsonEquationOfState(matTarget, etamin, etamax, units)
    eosTarget.minimumPressure = pmin
elif matTarget.lower() in shelpers.material_strings['m/aneos']:
    print "m/anoes" #TODO put in aneos
else:
    raise ValueError("invalid material selection for target")

# Select and construct impactor eos
if matImpactor.lower() in shelpers.material_strings['tillotson']:
    etamin, etamax = rhomin/rhoImpactor, rhomax/rhoImpactor
    pext, pmin, pmax = 0.0, 0.0, 1e200
    eosImpactor = TillotsonEquationOfState(matImpactor, etamin, etamax, units)
    eosImpactor.minimumPressure = pmin
elif matImpactor.lower() in shelpers.material_strings['m/aneos']:
    print "m/anoes" # TODO put in aneos
else:
    raise ValueError("invalid material selection for impactor")

# Verify valid EOSs
if eosTarget is None or eosImpactor is None:
    raise ValueError("target and/or impactor eos construction failed")
if not (eosTarget.valid() and eosImpactor.valid()):
    raise ValueError("target and/or impactor eos construction failed")

#-------------------------------------------------------------------------------
# NAV Restarts and output directories
# Here we create the output directories, and deal with restarted runs if any.
#-------------------------------------------------------------------------------
# Name directories and files.
jobDir = os.path.join(baseDir, 
                       'nxTarget=%i' % nxTarget,
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
# Here we create and populate a node list with initial conditions. In spheral, the
# construction order is as follows:
# 1. Create an empty node list with fields that match the size and type of problem.
# 2. Create a "generator" that will decide what values to give to all field variables
#    of node i. Normally we start with one of the simple, stock generators, and
#    modify the x,y,z,vx,vy,vz,rho,U values to suit our initial conditions.
# 3. Distribute, using the (nodeList, generator) pair, among ranks. The generator
#    will be used to fill values in the node list, and then discarded. 
#-------------------------------------------------------------------------------
# Create the node lists.
target   = makeFluidNodeList('target', eosTarget, 
                             nPerh = nPerh, 
                             xmin = -10.0*rTarget*Vector.one, # (probably unnecessary)
                             xmax =  10.0*rTarget*Vector.one, # (probably unnecessary)
                             hmin = hmin,
                             hmax = hmax,
                             rhoMin = rhomin,
                             rhoMax = rhomax,
                             )
impactor = makeFluidNodeList('impactor', eosImpactor, 
                             nPerh = nPerh, 
                             xmin = -10.0*rImpactor*Vector.one, # (probably unnecessary)
                             xmax =  10.0*rImpactor*Vector.one, # (probably unnecessary)
                             hmin = hmin,
                             hmax = hmax,
                             rhoMin = rhomin,
                             rhoMax = rhomax,
                             )
nodeSet = [target, impactor]

# Unless restarting, create the generators and set initial field values.
if restoreCycle is None:
    # Determine appropriate resolution for impactor.
    m_per_node_target = 1.0 * mTarget / (nxTarget**3)
    nxImp = max(2, int((mImpactor/m_per_node_target)**(1.0/3.0)))
    m_per_node_imp = 1.0 * mImpactor / (nxImp**3)
    print "Selected {} nodes across impactor.".format(nxImp)
    print "Target node mass = {}; Impactor node mass = {}".format(
                                                           m_per_node_target,
                                                           m_per_node_imp)

    # Start with the stock generators.
    targetGenerator   = GenerateNodeDistribution3d(nxTarget, nxTarget, nxTarget,
                                   rhoTarget,
                                   distributionType = 'lattice',
                                   xmin = (-rTarget, -rTarget, -rTarget),
                                   xmax = ( rTarget,  rTarget,  rTarget),
                                   rmin = 0.0,
                                   rmax = rTarget,
                                   nNodePerh = nPerh
                                   )
    impactorGenerator = GenerateNodeDistribution3d(nxImp, nxImp, nxImp,
                                   rhoImpactor,
                                   distributionType = 'lattice',
                                   xmin = (-rImpactor, -rImpactor, -rImpactor),
                                   xmax = ( rImpactor,  rImpactor,  rImpactor),
                                   rmin = 0.0,
                                   rmax = rImpactor,
                                   nNodePerh = nPerh
                                   )

    # Disturb the lattice symmetry to avoid artificial singularities.
    for k in range(targetGenerator.localNumNodes()):
        targetGenerator.x[k] *= 1.0 + random.uniform(-0.02, 0.02)
        targetGenerator.y[k] *= 1.0 + random.uniform(-0.02, 0.02)
        targetGenerator.z[k] *= 1.0 + random.uniform(-0.02, 0.02)
        pass
    for k in range(impactorGenerator.localNumNodes()):
        impactorGenerator.x[k] *= 1.0 + random.uniform(-0.02, 0.02)
        impactorGenerator.y[k] *= 1.0 + random.uniform(-0.02, 0.02)
        impactorGenerator.z[k] *= 1.0 + random.uniform(-0.02, 0.02)
        pass

    # Place the impactor at the point of impact. It is coming from the 
    # positive x direction in the xy plane.
    displace = Vector((rTarget+rImpactor)*cos(pi/180.0*angleImpact),
                      (rTarget+rImpactor)*sin(pi/180.0*angleImpact),
                      0.0)
    for k in range(impactorGenerator.localNumNodes()):
        impactorGenerator.x[k] += displace.x
        impactorGenerator.y[k] += displace.y
        impactorGenerator.z[k] += displace.z
                                                       
    # Fill node lists using generators and distribute to ranks.
    print "Starting node distribution..."
    distributeNodes3d((target, targetGenerator),
                      (impactor, impactorGenerator))
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

#-------------------------------------------------------------------------------
# NAV Spheral's simulation structure
# Here we construct the objects that compose spheral's simulation hierarchy.
# These are:
#  * One or more physics packages (hydro, gravity, strength, damage)
#  * A time integrator of some flavor (usually a Runge-Kutta 2)
#  * The simulation controller
#-------------------------------------------------------------------------------
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
#  * output() - a generic access routine, usually a pickle of node list or some
#               calculated value of interest [cycle or time based]
#-------------------------------------------------------------------------------
def mOutput(stepsSoFar,timeNow,dt):
    """Save node list to flat ascii file."""
    mFileName="{0}-{1:04d}-{2:g}.{3}".format(
              jobName, stepsSoFar, timeNow, 'fnl')
    shelpers.pflatten_node_list(target, outDir + '/' + mFileName)
    pass
if not outCycle is None:
    control.appendPeriodicWork(mOutput,outCycle)
if not outTime is None:
    control.appendPeriodicTimeWork(mOutput,outTime)

#-------------------------------------------------------------------------------
# NAV Launch simulation
# The simulation can be run for a specified number of steps, or a specified time
# in seconds.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
    control.dropRestartFile()
else:
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

#-------------------------------------------------------------------------------
# NAV Final thoughts
# Here we may print a message if desired, or do any final action.
#-------------------------------------------------------------------------------
print "\n", jobName.upper(), "completed."
