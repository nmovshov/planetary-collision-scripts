#! /auto/proj/nmovshov_hindmost/collisions/SPHERAL/bin/python
#-------------------------------------------------------------------------------
# Simulation of one sphere collapsing under self gravity.
# This script can serve as a template for equilibrating fluid planets using
# various equations of state.
# Copy this file to a separate folder and modify with appropriate initial
# conditions. Make sure the first line includes the full path to SPHERAL's
# python. Then run with mpirun.
#-------------------------------------------------------------------------------
from math import *
import sys, mpi
from findLastRestart import *
from VoronoiDistributeNodes import distributeNodes3d
from NodeHistory import NodeHistory

#-------------------------------------------------------------------------------
# NAV Identify job name here
#-------------------------------------------------------------------------------
jobName = "gravityCollapse"
jobDesc = "a uniform (initially) planet collapsing under self gravity"
print jobDesc

#-------------------------------------------------------------------------------
# NAV Setup
# Generic, mutable problem parameters, in MKS units.
# Modify these to suit your specific problem.
#-------------------------------------------------------------------------------

# Experiment geometry
rPlanet = 0.5             # (earth radii) initial radius of planet
mPlanet = 0.2             # (earth masses) initial mass of planet
matPlanet = "basalt"      # granite, basalt, nylon, pure ice, water
mPlanet *= 5.972e24
rPlanet *= 6371.0e3

# Gravity parameters
softLength = 1.0e-5       # (fraction of planet radius) softening length
opening = 1.0             # (dimensionless) opening parameter for gravity tree walk
fdt = 0.1                 # (dimensionless) gravity timestep multiplier
softLength *= rPlanet

# Node seeding parameters ("resolution")
nxTarget = 20             # Number of nodes across the diameter of the target
nPerh = 1.51              # Nominal number of nodes per smoothing scale

# Times, simulation control, and output
steps = None              # None or advance a number of steps rather than to a time
goalTime = 5000           # Time to advance to (sec)
dt = 100                  # Initial guess for time step (sec)
dtMin = 0.1               # Minimum allowed time step (sec)
dtMax = 10.0              # Maximum allowed time step (sec)
vizTime = 180             # Time frequency for dropping viz files (sec)
vizCycle = 800            # Cycle frequency for dropping viz files
baseDir = jobName         # Base name for directory to store output in

#-------------------------------------------------------------------------------
# NAV Additional global paremeters that are rarely changed
#-------------------------------------------------------------------------------

# More simulation parameters
dtGrowth = 2.0            # Maximum growth factor for time step in a cycle (dimensionless)
verbosedt = True          # Verbose reporting of the time step criteria per cycle
maxSteps = None           # Maximum allowed steps for simulation advance
statsStep = 10            # Frequency for sampling conservation statistics and such
redistributeStep = 1000   # Frequency to load balance problem from scratch
restartStep = 500         # Frequency to drop restart files
restoreCycle = None       # If restarting, cycle to start from (if None, latest available restart cycle is selected)

# Artificial viscosity (and other numerical crap).
HydroConstructor = SolidASPHHydro
Qconstructor = MonaghanGingoldViscosity
Cl = 1.0
Cq = 1.0
Qlimiter = False
balsaraCorrection = False
epsilon2 = 1e-2
negligibleSoundSpeed = 1e-1 #TODO make depend on physics
csMultiplier = 1e-4
hmin = 1.0e-3*rImpactor
hmax = 1.0e-1*rTarget
hminratio = 0.1
limitIdealH = False
cfl = 0.5
useVelocityMagnitudeForDt = False
XSPH = True
epsilonTensile = 0.3
nTensile = 4

# Hydro parameters.
HEvolution = IdealH                 # Algorithm for updating the H (smoothing scale) tensor
densityUpdate = IntegrateDensity    # Algorithm for updating mass density
compatibleEnergyEvolution = True    # Energy update choice (compatibleEnergyEvolution results in machine precision energy conservation)
rigorousBoundaries = False          # Do we re-compute ghost nodes during a timestep (more expensive if true)

#-------------------------------------------------------------------------------
# NAV Build EOS object
#-------------------------------------------------------------------------------
units = PhysicalConstants(1.0,   # Unit length in meters
                          1.0,   # Unit mass in kg
                          1.0)   # Unit time in seconds
etamin, etamax = 0.01, 100.0     # bounds of rho/rho0
eosPlanet = TillotsonEquationOfState(matTarget,etamin,etamax,units)

#-------------------------------------------------------------------------------
# NAV Here we compute some derived problem parameters
#-------------------------------------------------------------------------------
pass

#-------------------------------------------------------------------------------
# NAV Here we determine if, and deal with, restarted runs
#-------------------------------------------------------------------------------
# Restart and output files.
dataDir = os.path.join(baseDir, 
                       "nxTarget=%i" % nxTarget,
                       "angle_impact=%3.1f" % angle_impact,
                       "v_impact=%g_mpersec" % vImpact,
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
# NAV Here we construct a node list based on our problem's geometry and IC
#-------------------------------------------------------------------------------
# Create the NodeLists.
etaMin = 0.2
etaMax = 5.0
#TODO: makeFLuidNodeList
nodeSet = [planet]

#-------------------------------------------------------------------------------
# NAV PREPROCESS Here we set node properties (positions, velocites, etc.)
#-------------------------------------------------------------------------------
if restoreCycle is None:
    print "Generating node distribution."
    from GenerateNodeDistribution3d import GenerateNodeDistribution3d

    targetGenerator = GenerateNodeDistribution3d(nxTarget, nxTarget, nxTarget,
                                                 rhoTarget,
                                                 distributionType = "lattice",
                                                 xmin = (-rTarget, -rTarget, -rTarget),
                                                 xmax = ( rTarget,  rTarget,  rTarget),
                                                 rmax = rTarget,
                                                 nNodePerh = nPerh)
    impactorGenerator = GenerateNodeDistribution3d(nxImpactor, nxImpactor, nxImpactor,
                                                 rhoImpactor,
                                                 distributionType = "lattice",
                                                 xmin = (-rImpactor, -rImpactor, -rImpactor),
                                                 xmax = ( rImpactor,  rImpactor,  rImpactor),
                                                 rmax = rImpactor,
                                                 nNodePerh = nPerh)

    # The above logic generates node positions centered on (0,0,0), but we need
    # to displace the impactor so it is just touching the surface of the target
    # at the requested angle.  We'll have it coming in from the positive x direction
    # in the xy plane.
    disp = Vector((rImpactor + rTarget)*cos(pi/180.0*angle_impact),
                  (rImpactor + rTarget)*sin(pi/180.0*angle_impact),
                  0.0)
    for i in xrange(impactorGenerator.localNumNodes()):
        impactorGenerator.x[i] += disp.x
        impactorGenerator.y[i] += disp.y
        impactorGenerator.z[i] += disp.z

    print "Starting node distribution..."
    distributeNodes3d((target, targetGenerator),
                      (impactor,  impactorGenerator))

    nGlobalNodes = 0
    for n in nodeSet:
        print "Generator info for %s" % n.name
        print "   Minimum number of nodes per domain : ", mpi.allreduce(n.numInternalNodes, mpi.MIN)
        print "   Maximum number of nodes per domain : ", mpi.allreduce(n.numInternalNodes, mpi.MAX)
        print "               Global number of nodes : ", mpi.allreduce(n.numInternalNodes, mpi.SUM)
        nGlobalNodes += mpi.allreduce(n.numInternalNodes, mpi.SUM)
    del n
    print "Total number of (internal) nodes in simulation: ", nGlobalNodes
    print "Ratio of impactor/target node mass : ", impactor.mass().max()/target.mass().max()
    
    # Intialize the impactor velocity.
    vel = impactor.velocity()
    for i in xrange(impactor.numInternalNodes):
        vel[i].x = -vImpact

# Construct a DataBase to hold our node lists.
db = DataBase()
for n in nodeSet:
    db.appendNodeList(n)
del n

#-------------------------------------------------------------------------------
# NAV Here we create the various objects needed by spheral
#-------------------------------------------------------------------------------

# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
WT = TableKernel(BSplineKernel(), 1000)
WTPi = WT

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

# Construct a damage model.
if useDamage:
    damageModelTarget = DamageModelConstructor(target,
                                               kWeibull = kWeibull,
                                               mWeibull = mWeibull,
                                               kernel = WT,
                                               seed = randomSeed,
                                               volume = 0.0,  # forces internal computation.
                                               volumeStretchFactor = 1.0,
                                               strainAlgorithm = strainType,
                                               effectiveDamageAlgorithm = damageType,
                                               useDamageGradient = useDamageGradient,
                                               flawAlgorithm = effectiveFlawAlgorithm,
                                               criticalDamageThreshold = criticalDamageThreshold)

# Build some history objects to follow the time evolution of stuff.
    strainHistory = AverageStrain(damageModelTarget,
                                  os.path.join(dataDir, "strainhistory.txt"))

# Construct a time integrator.
integrator = CheapSynchronousRK2Integrator(db)
integrator.appendPhysicsPackage(hydro)
if useDamage:
    integrator.appendPhysicsPackage(damageModelTarget)
integrator.lastDt = dt
integrator.verbose = verbosedt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
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

# Register some of our diagnostics to be fired during advancement.
if useDamage:
    control.appendPeriodicWork(strainHistory.sample, strainFrequency)

    # If we restarted flush the history files to catch up with the current state.
    strainHistory.flushHistory()

#-------------------------------------------------------------------------------
# NAV MIDPROCESS Here we set register optional work to be done mid-run
#-------------------------------------------------------------------------------
def midprocess(stepsSoFar,timeNow,dt):
	pass
frequency=4000
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
    control.step() # One more step to ensure we get the final viz dump.

#-------------------------------------------------------------------------------
# NAV Here we do any post processing
#-------------------------------------------------------------------------------
#from IPython import embed
#embed() # uncomment to start an interactive session when the run completes
