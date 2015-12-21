#! /proj/nmovshov_hindmost/collisions/SPHERAL/bin/python
#-------------------------------------------------------------------------------
# Run two solid, strength-regime spheres into each other.
# 
# This script launches a basic strength regime collision. Two spherical, solid
# objects of arbitrary size and constant density collide with some specified
# velocity and impact parameter. In spheral terms, the target and impactor are
# of type SolidNodeList3d. The only physics package attached to the integrator
# is the hydro package. A damage model (GradyKippTensorDamageBenzAsphaug3d) is
# attached to the target node list.
#
# This script can be used to model collision of small, homogeneous bodies. 
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
jobName = 'solidcollision'
jobDesc = "Simple collision of solid, homogeneous spheres."
print '\n', jobName.upper(), '-', jobDesc.upper()

# Target parameters
rTarget = 30e3               # Target radius (m)
rhoTarget = 2700.0           # Target approximate bulk density (kg/m^3)
matTarget = 'basalt'         # Target material (see <pcs>/MATERIALS.md for options)
mTarget = 4.0/3.0*pi*rhoTarget*rTarget**3

# Impactor parameters
rImpactor = 30e3             # Impactor radius (m)
rhoImpactor = 2700.0         # Impactor approximate bulk density (kg/m^3)
matImpactor = 'basalt'       # Impactor material (see <pcs>/MATERIALS.md for options)
mImpactor = 4.0/3.0*pi*rhoImpactor*rImpactor**3

# Collision parameters
vImpact = 3000               # Impact velocity (m/s)
angleImpact = 2              # Impact angle to normal (degrees)
crossTime = 2*rTarget/vImpact

# Strength and damage parameters (TODO: bake into material choice!)
muTarget   = 2.27e10         # Shear modulus (Pa)
Y0Target   = 3.5e9           # Plastic yield stress (Pa)
muImpactor = 7.3e5           # Shear modulus (Pa)
Y0Impactor = 1.0e5           # Plastic yield stress (Pa)
kWeibull   = 8.0e43
mWeibull   = 9.5

# Times, simulation control, and output
nxTarget = 40                # Nodes across diameter of target (run "resolution")
steps = None                 # None or number of steps to advance (overrides time)
goalTime = 2*crossTime       # Time to advance to (sec)
dtInit = 0.02                # Initial guess for time step (sec)
vizTime = 0.2*goalTime       # Time frequency for dropping viz files (sec)
vizCycle = None              # Cycle frequency for dropping viz files
outTime = vizTime            # Time between running output routine (sec)
outCycle = None              # Cycles between running output routine

# Node list parameters
nPerh = 2.01                 # Nominal number of nodes per smoothing scale
hmin = 0.001                 # Minimum smoothing length (fraction of nominal)
hmax = 2.0                   # Maximum smoothing length (fraction of nominal)
rhomax = 1e+1*rhoTarget      # Upper bound on node density (kg/m^3)
generator_type = 'hcp'       # Node generator to use. 'hcp'|'old'|'shells'
hmin *= nPerh*2*rTarget/nxTarget
hmax *= nPerh*2*rTarget/nxTarget
rhomin = mTarget/nxTarget**3/hmax**3
universeEdge = 30*rTarget

# More simulation parameters
dtGrowth = 2.0               # Maximum growth factor for time step per cycle 
dtMin = 0                    # Minimum allowed time step (sec)
dtMax = 0.1*goalTime         # Maximum allowed time step (sec)
verbosedt = False            # Verbose reporting of the time step criteria 
maxSteps = 10000             # Maximum allowed steps for simulation advance
statsStep = None             # Frequency for sampling conservation statistics etc.
redistributeStep = 80000     # Frequency to load balance problem from scratch
restartStep = 200            # Frequency to drop restart files
restoreCycle = None          # If None, latest available restart cycle is selected
baseDir = jobName            # Base name for directory to store output in

# Damage and strength modeling (TODO: understand this)
DamageModelConstructor = GradyKippTensorDamageBenzAsphaug
randomSeed = 78987265
strainType = PseudoPlasticStrain
damageType = Copy
useDamageGradient = True
minFlawsPerNode = 1
effectiveFlawAlgorithm = SampledFlaws
criticalDamageThreshold = 0.5

# Cooldown mechanism (normally disabled, with cooldownFrequency = None)
cooldownMethod = 'dashpot'   # 'dashpot' or 'stomp' 
cooldownPower = 0.1          # Dimensionless cooldown "strength" >=0
cooldownFrequency = None     # Cycles between application (use 1 with dashpot)
                             # * With 'stomp' method, 0<=power<=1

# Node Culling mechanism (normally disabled, with cullingFrequency = None)
cullingFrequency = None      # Cycles between application

#-------------------------------------------------------------------------------
# NAV Assertions
# This is a good place for a quick abort if some bad parameter choices are going
# to cause trouble later, in confusing ways. We assume that spheral constructors
# use their own assertions, so here we can validate just our own stuff. Another
# valid option would be to simply not worry about it, and let exceptions happen.
#-------------------------------------------------------------------------------
assert rTarget >= rImpactor
assert vImpact >= 0.
assert 0 <= angleImpact < 90, "give impact angle in first quadrant (in degrees)"
assert (outTime is None) or (outCycle is None),\
        "output on both time and cycle is confusing"
assert generator_type in ['hcp', 'shells', 'old']
if cooldownFrequency is not None:
    sys.stderr.write("\033[1;31m")
    print "WARNING - damping is enabled, is this intentional?"
    sys.stderr.write("\033[0m")
    assert 0 <= cooldownPower, "cool DOWN not up"
    if cooldownMethod is 'stomp':
        assert 0 <= cooldownPower <= 1.0, "stomp fraction is 0-1"
    assert type(cooldownFrequency) is int and cooldownFrequency > 0, "very funny"
    assert cooldownMethod in ['dashpot','stomp'], "unknown cooldown method"
    assert (cooldownFrequency == 1) or (not(cooldownMethod is 'dashpot')),\
            "dashpot cooling method requires frequency=1"
if cullingFrequency is not None:
    sys.stderr.write("\033[1;31m")
    print "WARNING - node culling is enabled, is this intentional?"
    sys.stderr.write("\033[0m")
    assert type(cullingFrequency) is int and cullingFrequency > 0

#-------------------------------------------------------------------------------
# NAV Spheral hydro solver options
# These options for spheral's hydro mechanism are normally left alone.
#-------------------------------------------------------------------------------
HydroConstructor = SolidASPHHydro
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
# Here we construct equation-of-state objects, one per node list. In this case,
# one each for the target and impactor. The choice of eos is determined by the 
# material string (see <pcs>MATERIALS.md for available options).
#-------------------------------------------------------------------------------
eosTarget, eosImpactor = None, None

# Most eos constructors need to know about units. We usually use MKS.
units = PhysicalConstants(1.0, # unit length in meters
                          1.0, # unit mass in kilograms
                          1.0) # unit time in seconds

# Construct and verify target eos
eosTarget = shelpers.construct_eos_for_material(matTarget,units)
assert eosTarget is not None
assert eosTarget.valid()

# Construct and verify impactor eos
eosImpactor = shelpers.construct_eos_for_material(matImpactor,units)
assert eosImpactor is not None
assert eosImpactor.valid()

#-------------------------------------------------------------------------------
# NAV Prepare strength and damage models for target and impactor
#-------------------------------------------------------------------------------
strengthTarget =   ConstantStrength(muTarget,
                                    Y0Target)
strengthImpactor = ConstantStrength(muImpactor,
                                    Y0Impactor)

# Optionally, provide non-default values to the following
eosTarget.etamin_solid = 0.94 # default is 0.94
eosTarget.minimumPressure = 0.0 # default is -1e+200
eosImpactor.etamin_solid = 0.94 # default is 0.94
eosImpactor.minimumPressure = 0.0 # default is -1e+200

#-------------------------------------------------------------------------------
# NAV Restarts and output directories
# Here we create the output directories, and deal with restarted runs if any.
#-------------------------------------------------------------------------------
# Name directories and files.
jobDir = os.path.join(baseDir, 
                       'rTarget=%0.2g' % rTarget,
                       'eosTarget=%d' % eosTarget.uid,
                       'rImpactor=%0.2g' % rImpactor,
                       'eosImpactor=%d' % eosImpactor.uid,
                       'angle=%g' % angleImpact,
                       'vImpact=%0.3g' % vImpact,
                       'nxTarget=%i' % nxTarget,
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
target   = makeSolidNodeList('target', eosTarget, strengthTarget,
                             nPerh = nPerh, 
                             xmin = -universeEdge*Vector.one, # (probably unnecessary)
                             xmax =  universeEdge*Vector.one, # (probably unnecessary)
                             hmin = hmin,
                             hmax = hmax,
                             rhoMin = rhomin,
                             rhoMax = rhomax,
                             hminratio = hminratio,
                             topGridCellSize = universeEdge,
                             )
target.eos_id = eosTarget.uid

impactor = makeSolidNodeList('impactor', eosImpactor, strengthImpactor,
                             nPerh = nPerh, 
                             xmin = -universeEdge*Vector.one, # (probably unnecessary)
                             xmax =  universeEdge*Vector.one, # (probably unnecessary)
                             hmin = hmin,
                             hmax = hmax,
                             rhoMin = rhomin,
                             rhoMax = rhomax,
                             hminratio = hminratio,
                             topGridCellSize = universeEdge,
                             )
impactor.eos_id = eosImpactor.uid

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

    # Create a basic, usually constant density generator.
    if generator_type == 'old':
        targetGenerator   = GenerateNodeDistribution3d(nxTarget, nxTarget, nxTarget,
                              rhoTarget,
                              distributionType = 'lattice',
                              xmin = (-rTarget, -rTarget, -rTarget),
                              xmax = ( rTarget,  rTarget,  rTarget),
                              rmin = 0.0,
                              rmax = rTarget,
                              nNodePerh = nPerh)
        impactorGenerator = GenerateNodeDistribution3d(nxImp, nxImp, nxImp,
                              rhoImpactor,
                              distributionType = 'lattice',
                              xmin = (-rImpactor, -rImpactor, -rImpactor),
                              xmax = ( rImpactor,  rImpactor,  rImpactor),
                              rmin = 0.0,
                              rmax = rImpactor,
                              nNodePerh = nPerh)
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
        pass
    elif generator_type == 'hcp':
        targetGenerator   = PlanetNodeGenerators.HexagonalClosePacking(
                              nx = nxTarget,
                              rho = rhoTarget,
                              scale = 2*rTarget,
                              rMin = 0.0,
                              rMax = rTarget,
                              nNodePerh = nPerh)
        impactorGenerator = PlanetNodeGenerators.HexagonalClosePacking(
                              nx = nxImp,
                              rho = rhoImpactor,
                              scale = 2*rImpactor,
                              rMin = 0.0,
                              rMax = rImpactor,
                              nNodePerh = nPerh)
        pass
    elif generator_type == 'shells':
        targetGenerator   = PlanetNodeGenerators.EqualSpacingSphericalShells(
                              nLayers = nxTarget/2,
                              rho = rhoTarget,
                              rMin = 0.0,
                              rMax = rTarget,
                              nNodePerh = nPerh)
        impactorGenerator = PlanetNodeGenerators.EqualSpacingSphericalShells(
                              nLayers = nxImp/2,
                              rho = rhoImpactor,
                              rMin = 0.0,
                              rMax = rImpactor,
                              nNodePerh = nPerh)
        pass
    else:
        print "unknown generator type"
        sys.exit(1)
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
        pass
                                                       
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
    WMR = (max(impactor.mass().max(), target.mass().max())/
           min(impactor.mass().min(), target.mass().min()))
    if WMR < 1.5:
        sys.stderr.write("\033[32m")
    else:
        sys.stderr.write("\033[31m")
    print "Worst node mass ratio: {}".format(WMR)
    sys.stderr.write("\033[0m")
    assert WMR < 3, "Severe node mass imbalance: comment out assertion to ignore."

    # Launch the impactor
    vel = impactor.velocity()
    for k in range(impactor.numInternalNodes):
        vel[k].x = -vImpact
        pass

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

# Construct a damage model. (TODO: understand this and modify to use material tag)
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

# Create the time integrator and attach the physics packages to it.
integrator = CheapSynchronousRK2Integrator(db)
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
#  * cullnodes() - remove errant nodes from simulation (use with care)
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

def cullnodes(stepsSoFar,timeNow,dt):
    for nl in nodeSet:
        v = nl.velocity()
        r = nl.positions()
        u = nl.specificThermalEnergy()
        C = db.connectivityMap()
        bads = vector_of_int()
        for k in range(nl.numInternalNodes):
            if ( abs(v[k].magnitude()) > 10*vImpact   or
                 u[k] > 10*eosTarget.epsVapor         or
                 0 < C.numNeighborsForNode(nl,k) < 10 or
                 sqrt(r[k].x**2 + r[k].y**2 + r[k].z**2) > universeEdge
               ):
                bads.append(k)
                pass
            pass
        nl.deleteNodes(bads)
        if bads.size() > 0:
            sys.stderr.write("\033[31m")
            msg = "WARNING - cycle {} deleted {} nodes from {}\n".format(stepsSoFar,bads.size(),nl.name)
            sys.stderr.write(msg)
            sys.stderr.write("\033[0m")
        pass
    pass
control.appendPeriodicWork(cullnodes,cullingFrequency)

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

#-------------------------------------------------------------------------------
# NAV Final thoughts
# Here we may print a message if desired, or do any final action.
#-------------------------------------------------------------------------------
print "\n", jobName.upper(), "completed."
