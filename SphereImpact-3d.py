#-------------------------------------------------------------------------------
# Simulation of the sphere colliding experiments described in
#
# Nakamura, A. and A. Fujiwara, 1991, Velocity Distribution of fragments Formed 
# in a Simulated Collisional Disruption. Icarus, 92, 132-146.
#-------------------------------------------------------------------------------
from math import *
import sys, mpi
from SolidSpheral3d import *
from findLastRestart import *

from VoronoiDistributeNodes import distributeNodes3d

from NodeHistory import NodeHistory
from AverageStrain import AverageStrain

#-------------------------------------------------------------------------------
# Identify ourselves!
#-------------------------------------------------------------------------------
print "Nakamura & Fujiwara colliding spheres experiment."

#-------------------------------------------------------------------------------
# Generic problem parameters.
# All (centimeter, gram, microsec) units.
#
# Note that variables passed in the "commandLine" function can be overridden on
# the command line, i.e.,
#       "python SphereImpact-3d.py --nxTarget 100 --mImpactor 0.5"
#-------------------------------------------------------------------------------
commandLine(
    # Experiment geometry
    rTarget = 3.0,                   # (cm) Radius of target
    mImpactor = 0.2,                 # (g) mass of impactor
    vImpactor = 0.32,                # (cm/usec) initial velocity of impactor
    angle_impactor = 30.0,           # Impact angle to normal (degrees)
    
    # Node seeding parameters.
    nxTarget = 40,                   # Number of nodes across the diameter of the target
    nPerh = 1.51,                    # Nominal number of nodes per smoothing scale

    # Material properties
    # Given here as cm-gm-usec units.

    # Gruneisen parameters for Lucite (used as a stand-in for Nylon).
    rho0_Nylon = 1.185,               # (g/cm^3)
    C0_Nylon = 0.218,                # (cm/usec)
    S1_Nylon = 2.08,                 # (dimensionless)
    S2_Nylon = -1.124,               # (dimensionless)
    S3_Nylon = 0.0,                  # (dimensionless)
    gamma0_Nylon = 0.85,             # (dimensionless)
    b_Nylon = 0.0,                   # (dimensionless)
    atomicWeight_Nylon = 226.32,     # atomic weight
    mu_Nylon = 7.3e-4,               # Shear modulus
    Y0_Nylon = 1.0e-4,               # Plastic yield stress
    etaMin_Nylon = 0.5,              # min rho/rho0 (dimensionless)
    etaMax_Nylon = 5.0,              # max rho/rho0 (dimensionless)

    # Basalt 
    rho0_Basalt = 2.70,              # (g/cm^3)
    a_Basalt = 0.5,                  # Tillotson, dimensionless
    b_Basalt = 1.5,                  # Tillotson, dimensionless
    A_Basalt = 0.267,                # Tillotson  (pressure = g/(cm usec^2))
    B_Basalt = 0.267,                # Tillotson, (pressure = g/(cm usec^2))
    alpha_Basalt = 5.0,              # Tillotson, dimensionless
    beta_Basalt = 5.0,               # Tillotson, dimensionless
    eps0_Basalt = 4.87,              # Tillotson, ( (cm/usec)^2 )
    epsLiquid_Basalt = 4.72e-2,      # Tillotson, ( (cm/usec)^2 )
    epsVapor_Basalt = 1.82e-1,       # Tillotson, ( (cm/usec)^2 )
    atomicWeight_Basalt = 60.08,     # Tillotson, atomic weight
    mu_Basalt = 2.27e-1,             # Shear modulus
    Y0_Basalt = 3.5e-2,              # Plastic yield stress
    kWeibull_Basalt = 5.0e28,        # Damage model, cm^-3
    mWeibull_Basalt = 8.5,           # Damage model, dimensionless
    etaMin_Basalt = 0.2,             # min rho/rho0 (dimensionless)
    etaMax_Basalt = 5.0,             # max rho/rho0 (dimensionless)

    # Switches to turn off physics.
    useStrength = True,
    useDamage = True,

    # Damage and strength modeling
    DamageModelConstructor = GradyKippTensorDamageBenzAsphaug,
    randomSeed = 78987265,
    strainType = PseudoPlasticStrain,
    damageType = Copy,
    useDamageGradient = True,
    minFlawsPerNode = 1,
    effectiveFlawAlgorithm = SampledFlaws,
    criticalDamageThreshold = 0.5,

    # Artificial viscosity (and other numerical crap).
    HydroConstructor = SolidASPHHydro,
    Qconstructor = MonaghanGingoldViscosity,
    Cl = 1.0,
    Cq = 1.0,
    Qlimiter = False,
    balsaraCorrection = False,
    epsilon2 = 1e-2,
    negligibleSoundSpeed = 1e-5,
    csMultiplier = 1e-4,
    hmin = 1.0e-2,
    hmax = 1e3,
    hminratio = 0.1,
    limitIdealH = False,
    cfl = 0.5,
    useVelocityMagnitudeForDt = False,
    XSPH = True,
    epsilonTensile = 0.3,
    nTensile = 4,

    # Hydro parameters.
    HEvolution = IdealH,                 # Algorithm for updating the H (smoothing scale) tensor
    densityUpdate = IntegrateDensity,    # Algorithm for updating mass density
    compatibleEnergyEvolution = True,    # Energy update choice (compatibleEnergyEvolution results in machine precision energy conservation)
    rigorousBoundaries = False,          # Do we re-compute ghost nodes at during a timestep (more expensive if true)

    # Times, and simulation control.
    steps = None,              # Optionally advance a number of steps rather than to a time
    goalTime = 100.0,           # Time to advance to (usec)
    dt = 1.0e-3,               # Initial guess for time step (usec)
    dtMin = 1e-5,              # Minimum allowed time step (usec)
    dtMax = 10.0,              # Maximum allowed time step (usec)
    dtGrowth = 2.0,            # Maximum growth factor for time step in a cycle (dimensionless)
    verbosedt = False,         # Verbose reporting of the time step criteria per cycle
    maxSteps = None,           # Maximum allowed steps for simulation advance
    statsStep = 10,            # Frequency for sampling conservation statistics and such
    redistributeStep = 1000,   # Frequency to load balance problem from scratch
    restartStep = 500,          # Frequency to drop restart files
    restoreCycle = None,       # If restarting, cycle to start from (if None, latest available restart cycle is selected)

    # Output
    vizTime = 10.,                     # Time frequency for dropping viz files (usec)
    vizCycle = 5000,                     # Cycle frequency for dropping viz files
    strainFrequency = 10,              # Cycle frequency for measuring strain in the target
    baseDir = "SphereImpact-output",   # Base name for directory to store output in
    )

#-------------------------------------------------------------------------------
# Get the variables straightened out
#-------------------------------------------------------------------------------
# Our unit set choice.
units = PhysicalConstants(0.01,    # Unit length in meters
                          0.001,   # Unit mass in kilograms
                          1.0e-6)  # Unit time in sec
# Determine the impactor resolution based on mass matching the target resolution.
m_per_point_Target = (2.0*rTarget/nxTarget)**3 * rho0_Basalt
rImpactor = (mImpactor/(4.0/3.0*pi*rho0_Nylon))**(1.0/3.0)
nxImpactor = max(2, int((8.0*rImpactor**3 * rho0_Nylon/m_per_point_Target)**(1.0/3.0)*sqrt(1.0) + 0.5))
print "Selected %i nodes across diameter of impactor of radius %g for target point mass of %g." % (nxImpactor, rImpactor, m_per_point_Target)

# You have to use strength if you're applying damage.
if useDamage:
    assert useStrength

#------------------------------------------------------------------------------
# Restart and output files.
#------------------------------------------------------------------------------
dataDir = os.path.join(baseDir, 
                       "nxTarget=%i" % nxTarget,
                       "angle_impact=%3.1f" % angle_impactor,
                       "v_impact=%g_kmpersec" % (vImpactor * 10.0),
                       )

restartDir = os.path.join(dataDir, "restarts", "proc-%04i" % mpi.rank)
vizDir = os.path.join(dataDir, "viz")
baseName = "SphereImpact-3d"
restartName = os.path.join(restartDir, baseName)

#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
# If we're restarting, find the set of most recent restart files.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    restoreCycle = findLastRestart(restartName)

#-------------------------------------------------------------------------------
# Basalt material properties.
#-------------------------------------------------------------------------------
eosBasalt = TillotsonEquationOfState(rho0_Basalt,        # Ref density
                                     etaMin_Basalt,      # (dimensionless)
                                     etaMax_Basalt,      # (dimensionless)
                                     a_Basalt,           # (dimensionless)
                                     b_Basalt,           # (dimensionless)
                                     A_Basalt,           # (pressure units)
                                     B_Basalt,           # (pressure units)
                                     alpha_Basalt,       # (dimensionless)
                                     beta_Basalt,        # (dimensionless)
                                     eps0_Basalt,        # (energy/mass)
                                     epsLiquid_Basalt,   # (energy/mass)
                                     epsVapor_Basalt,    # (energy/mass)
                                     atomicWeight_Basalt,# (dimensionless)
                                     units,
                                     )

if useStrength:
    strengthBasalt = ConstantStrength(mu_Basalt,   # shear modulus (Pa)
                                      Y0_Basalt)   # plastic yield stress (Pa)
else:
    strengthBasalt = NullStrength()

#-------------------------------------------------------------------------------
# Nylon material properties (though we are actually using Lucite)
#-------------------------------------------------------------------------------
eosNylon = GruneisenEquationOfState(rho0_Nylon,        # Ref density
                                    etaMin_Nylon,      # (dimensionless)
                                    etaMax_Nylon,      # (dimensionless)
                                    C0_Nylon,          # (velocity units)
                                    S1_Nylon,          # (dimensionless)
                                    S2_Nylon,          # (dimensionless)
                                    S3_Nylon,          # (dimensionless)
                                    gamma0_Nylon,      # (dimensionless)
                                    b_Nylon,           # (dimensionless)
                                    atomicWeight_Nylon,# (dimensionless)
                                    units,
                                    )

if useStrength:
    strengthNylon = ConstantStrength(mu_Nylon,   # shear modulus (Pa)
                                     Y0_Nylon)   # plastic yield stress (Pa)
else:
    strengthNylon = NullStrength()

#-------------------------------------------------------------------------------
# Create the NodeLists.
#-------------------------------------------------------------------------------
target = makeSolidNodeList("Target", eosBasalt, strengthBasalt,
                           nPerh = nPerh,
                           hmin = hmin,
                           hmax = hmax,
                           rhoMin = etaMin_Basalt*rho0_Basalt,
                           rhoMax = etaMax_Basalt*rho0_Basalt,
                           xmin = -1e2 * Vector.one,
                           xmax =  1e2 * Vector.one)
impactor = makeSolidNodeList("Impactor", eosNylon, strengthNylon,
                             nPerh = nPerh,
                             hmin = hmin,
                             hmax = hmax,
                             rhoMin = etaMin_Nylon*rho0_Nylon,
                             rhoMax = etaMax_Nylon*rho0_Nylon,
                             xmin = -1e2 * Vector.one,
                             xmax =  1e2 * Vector.one)
nodeSet = [target, impactor]

#-------------------------------------------------------------------------------
# Set node properties (positions, velocites, etc.)
#-------------------------------------------------------------------------------
if restoreCycle is None:
    print "Generating node distribution."
    from GenerateNodeDistribution3d import GenerateNodeDistribution3d

    targetGenerator = GenerateNodeDistribution3d(nxTarget, nxTarget, nxTarget,
                                                 rho0_Basalt,
                                                 distributionType = "lattice",
                                                 xmin = (-rTarget, -rTarget, -rTarget),
                                                 xmax = ( rTarget,  rTarget,  rTarget),
                                                 rmax = rTarget,
                                                 nNodePerh = nPerh)
    impactorGenerator = GenerateNodeDistribution3d(nxImpactor, nxImpactor, nxImpactor,
                                                 rho0_Nylon,
                                                 distributionType = "lattice",
                                                 xmin = (-rImpactor, -rImpactor, -rImpactor),
                                                 xmax = ( rImpactor,  rImpactor,  rImpactor),
                                                 rmax = rImpactor,
                                                 nNodePerh = nPerh)

    # The above logic generates node positions centered on (0,0,0), but we need
    # to displace the impactor so it is just touching the surface of the target
    # at the requested angle.  We'll have it coming in from the positive x direction
    # in the xy plane.
    disp = Vector((rImpactor + rTarget)*cos(pi/180.0*angle_impactor),
                  (rImpactor + rTarget)*sin(pi/180.0*angle_impactor),
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
        vel[i].x = -vImpactor

#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node lists.
#-------------------------------------------------------------------------------
db = DataBase()
for n in nodeSet:
    db.appendNodeList(n)
del n

#-------------------------------------------------------------------------------
# Create our interpolation kernels -- one for normal hydro interactions, and
# one for use with the artificial viscosity
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
WTPi = WT

#-------------------------------------------------------------------------------
# Construct the artificial viscosity.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
q.epsilon2 = epsilon2
q.negligibleSoundSpeed = negligibleSoundSpeed
q.csMultiplier = csMultiplier

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
# Construct a damage model.
#-------------------------------------------------------------------------------
# These terms to help hoke up the damage number so that we pretend the asteroid
# is actually a cylinder extruded in the z-direction by rOuter.
if useDamage:
    damageModelTarget = DamageModelConstructor(target,
                                               kWeibull = kWeibull_Basalt,
                                               mWeibull = mWeibull_Basalt,
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

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
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
# Register some of our diagnostics to be fired during advancement.
#-------------------------------------------------------------------------------
if useDamage:
    control.appendPeriodicWork(strainHistory.sample, strainFrequency)

    # If we restarted flush the history files to catch up with the current state.
    strainHistory.flushHistory()

#-------------------------------------------------------------------------------
# Advance to completion.
#-------------------------------------------------------------------------------
if not steps is None:
    control.step(steps)
    #raise ValueError, ("Completed %i steps." % steps)

else:
    control.advance(goalTime, maxSteps)
    control.dropRestartFile()
    control.step() # One more step to ensure we get the final viz dump.
