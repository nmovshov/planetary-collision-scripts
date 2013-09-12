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
rTarget = 200              # Target radius (m)
mTarget = 200              # Target mass (kg)
matTarget = 'basalt'       # Target material (see uss/MATERIALS.md for options)
rhoTarget = 3.0*mTarget/(4.0*pi*rTarget**3)

# Impactor parameters
rImpactor = 200            # Impactor radius (m)
mImpactor = 200            # Impactor mass (kg)
matImpactor = 'basalt'     # Impactor material (see uss/MATERIALS.md for options)
rhoImpactor = 3.0*mImpactor/(4.0*pi*rImpactor**3)

# Collision parameters
vImpact = 10               # Impact velocity (m/s)
angleImpact = 30           # Impact angle to normal (degrees)

# Times, simulation control, and output
nxTarget = 20              # Nodes across diameter of target (run "resolution")
steps = None               # None or advance a number of steps rather than to a time
goalTime = 6               # Time to advance to (sec)
dtInit = 0.02              # Initial guess for time step (sec)
vizTime = None             # Time frequency for dropping viz files (sec)
vizCycle = None            # Cycle frequency for dropping viz files
outTime = None             # Time between running output routine (sec)
outCycle = None            # Cycles between running output routine

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
    etamin, etamax = rhomin/rhoTarget, rhomin/rhoTarget
    eosTarget = TillotsonEquationOfState(matTarget, etamin, etamax, units)
elif matTarget.lower() in shelpers.material_strings['m/aneos']:
    print "m/anoes" #TODO put in aneos
else:
    raise ValueError("invalid material selection for target")

# Select and construct impactor eos
if matImpactor.lower() in shelpers.material_strings['tillotson']:
    etamin, etamax = rhomin/rhoImpactor, rhomax/rhoImpactor
    eosImpactor = TillotsonEquationOfState(matImpactor, etamin, etamax, units)
elif matImpactor.lower() in shelpers.material_strings['m/aneos']:
    print "m/anoes" # TODO put in aneos
else:
    raise ValueError("invalid material selection for impactor")

# Verify valid EOSs
if eosTarget is None or eosImpactor is None:
    raise ValueError("target and/or impactor eos construction failed")
if not (eosTarget.valid() and eosImpactor.valid()):
    raise ValueError("target and/or impactor eos construction failed")
