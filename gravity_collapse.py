#! /auto/proj/nmovshov_hindmost/collisions/SPHERAL/bin/python
#-------------------------------------------------------------------------------
# Quasi static collapse of a spherical planet.
# This script can serve as a template for equilibrating fluid planets using
# various equations of state.
# Copy this file to a separate folder and modify with appropriate initial
# conditions. Make sure the first line includes the full path to SPHERAL's
# python. Then run with mpirun.
#-------------------------------------------------------------------------------
from math import *
import sys, mpi
from Spheral3d import *
from findLastRestart import *
from VoronoiDistributeNodes import distributeNodes3d
from NodeHistory import NodeHistory

#-------------------------------------------------------------------------------
# NAV Identify job name here
#-------------------------------------------------------------------------------
jobName = "gravityCollapse"
jobDesc = "Quasi-static collapse of a uniform density (initially) planet."
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
rhoPlanet = 3*mPlanet/(4*pi*rPlanet**3)

# Gravity parameters
softLength = 1.0e-5       # (fraction of planet radius) softening length
opening = 1.0             # (dimensionless) opening parameter for gravity tree walk
fdt = 0.1                 # (dimensionless) gravity timestep multiplier
softLength *= rPlanet
G = MKS().G

# Node seeding parameters ("resolution")
nxPlanet = 20             # Number of nodes across the diameter of the target
nPerh = 1.51              # Nominal number of nodes per smoothing scale
hmin = 1.0e-3*rPlanet     # Lower bound on smoothing length
hmax = 1.0e-1*rPlanet     # Upper bound on smoothing length

# Times, simulation control, and output
steps = None              # None or advance a number of steps rather than to a time
goalTime = 5000.0         # Time to advance to (sec)
dt = goalTime/200         # Initial guess for time step (sec)
dtMin = 0.001*dt          # Minimum allowed time step (sec)
dtMax = 1000.0*dt         # Maximum allowed time step (sec)
vizTime = goalTime/20     # Time frequency for dropping viz files (sec)
vizCycle = 800            # Cycle frequency for dropping viz files
baseDir = jobName         # Base name for directory to store output in

# More simulation parameters
dtGrowth = 2.0            # Maximum growth factor for time step in a cycle (dimensionless)
verbosedt = False         # Verbose reporting of the time step criteria per cycle
maxSteps = None           # Maximum allowed steps for simulation advance
statsStep = 10            # Frequency for sampling conservation statistics and such
redistributeStep = 400    # Frequency to load balance problem from scratch
restartStep = 600         # Frequency to drop restart files
restoreCycle = None       # If None, latest available restart cycle is selected

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
negligibleSoundSpeed = 1e-4
csMultiplier = 1e-4
hminratio = 0.1
limitIdealH = False
cfl = 0.5
useVelocityMagnitudeForDt = False
XSPH = True
epsilonTensile = 0.3
nTensile = 4
HEvolution = IdealH
densityUpdate = IntegrateDensity
compatibleEnergyEvolution = True
rigorousBoundaries = False

#-------------------------------------------------------------------------------
# NAV Build EOS object
#-------------------------------------------------------------------------------
units = PhysicalConstants(1.0,   # Unit length in meters
                          1.0,   # Unit mass in kg
                          1.0)   # Unit time in seconds
etamin, etamax = 0.01, 100.0     # bounds of rho/rho0
#eosPlanet = TillotsonEquationOfState(matTarget,etamin,etamax,units)
eosPlanet = GammaLawGasMKS3d(gamma = 5.0/3.0, mu = 1.0)

#-------------------------------------------------------------------------------
# NAV Here we compute some derived problem parameters
#-------------------------------------------------------------------------------
pass

#-------------------------------------------------------------------------------
# NAV Here we determine if, and deal with, restarted runs
#-------------------------------------------------------------------------------
# Restart and output files.
dataDir = os.path.join(baseDir, 
                       "rPlanet=%g_m" % rPlanet,
                       "mPlanet=%g_kg" % mPlanet,
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
# NAV Here we construct a node list based on our problem's geometry and IC
#-------------------------------------------------------------------------------
# Create the NodeList.
planet = makeFluidNodeList("planet", eosPlanet, 
                           nPerh = nPerh, 
                           xmin = -10.0*rPlanet*Vector.one,
                           xmax =  10.0*rPlanet*Vector.one,
                           hmin = hmin,
                           hmax = hmax,
                           )
nodeSet = [planet]

#-------------------------------------------------------------------------------
# NAV PREPROCESS Here we set node properties (positions, velocites, etc.)
#-------------------------------------------------------------------------------
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

    # The above logic generates node positions centered on (0,0,0). Modify 
    # positions if necessary below.
    for k in range(planetGenerator.localNumNodes()):
        planetGenerator.x[k] += 0.0
        planetGenerator.y[k] += 0.0
        planetGenerator.z[k] += 0.0

    # Distribute nodes across ranks
    print "Starting node distribution..."
    distributeNodes3d((planet, planetGenerator))
    nGlobalNodes = 0
    for n in nodeSet:
        print "Generator info for %s" % n.name
        print "   Minimum number of nodes per domain : ", mpi.allreduce(n.numInternalNodes, mpi.MIN)
        print "   Maximum number of nodes per domain : ", mpi.allreduce(n.numInternalNodes, mpi.MAX)
        print "               Global number of nodes : ", mpi.allreduce(n.numInternalNodes, mpi.SUM)
        nGlobalNodes += mpi.allreduce(n.numInternalNodes, mpi.SUM)
    del n
    print "Total number of (internal) nodes in simulation: ", nGlobalNodes
    
    # Give initial velocity if desired.
    vel = planet.velocity()
    for k in range(planet.numInternalNodes):
        vel[k].x = 0.0

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
