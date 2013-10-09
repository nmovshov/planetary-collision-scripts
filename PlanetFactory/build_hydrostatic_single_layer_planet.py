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
jobName = 'onelayer'
jobDesc = "Hydrostatic equilibrium of a single-layer, fluid planet."
print '\n', jobName.upper(), '-', jobDesc.upper()

# Planet parameters
rPlanet = 500e3              # Initial guess for planet radius (m)
mPlanet = 5e22               # Total (and conserved) planet mass (kg)
matPlanet = 'basalt'         # Planet material (see <uss>/MATERIALS.md for options)
rhoPlanet = 3.0*mPlanet/(4.0*pi*rPlanet**3)
gravTime = 1/sqrt(MKS().G*rhoPlanet)

# Cooldown mechanism
cooldownMethod = 'dashpot'   # 'dashpot' or 'stomp' 
cooldownPower = 1.0          # Dimensionless cooldown "strength" >=0
cooldownFrequency = 1        # Cycles between application (use 1 with dashpot)
                             # * With 'stomp' method, 0<=power<=1

# Times, simulation control, and output
nxPlanet = 20                # Nodes across diameter of planet (run "resolution")
steps = None                 # None or advance a number of steps rather than to a time
goalTime = 2                 # Time to advance to (sec)
dtInit = 0.02                # Initial guess for time step (sec)
vizTime = 1                  # Time frequency for dropping viz files (sec)
vizCycle = None              # Cycle frequency for dropping viz files
outTime = None               # Time between running output routine (sec)
outCycle = None              # Cycles between running output routine

# Node list parameters
nPerh = 1.51                 # Nominal number of nodes per smoothing scale
hmin = int(1e-6*rPlanet)     # Lower bound on smoothing length
hmax = int(2e+0*rPlanet)     # Upper bound on smoothing length
rhomin = int(1e-6*rhoPlanet) # Lower bound on node density
rhomax = int(1e+1*rhoPlanet) # Upper bound on node density

# Gravity parameters
softLength = 1.0/nxPlanet    # Fraction of planet radius as softening length
opening = 1.0                # Dimensionless opening parameter for gravity tree walk
fdt = 0.1                    # Time step multiplier (dt=fdt*sqrt(softlength/a))
softLength *= rPlanet
G = MKS().G

# More simulation parameters
dtGrowth = 2.0               # Maximum growth factor for time step per cycle (dimensionless)
dtMin = 0                    # Minimum allowed time step (sec)
dtMax = int(0.1*goalTime)    # Maximum allowed time step (sec)
verbosedt = False            # Verbose reporting of the time step criteria per cycle
maxSteps = 800               # Maximum allowed steps for simulation advance
statsStep = None             # Frequency for sampling conservation statistics and such
redistributeStep = 2000      # Frequency to load balance problem from scratch
restartStep = 200            # Frequency to drop restart files
restoreCycle = None          # If None, latest available restart cycle is selected
baseDir = jobName            # Base name for directory to store output in
