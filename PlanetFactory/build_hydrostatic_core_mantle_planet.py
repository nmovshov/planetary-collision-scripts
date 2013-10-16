#! /proj/nmovshov_hindmost/collisions/SPHERAL/bin/python
#-------------------------------------------------------------------------------
# Set up a two-layer, fluid planet and run to hydrostatic equilibrium.
#
# This script serves as a template for equilibrating a fluid, two-layer  planet,
# specifying core and mantle radius and density.
# Copy this file to a separate folder and modify with choices of layer densities,
# radii, and materials. The parameter nxPlanet controls the run resolution. You
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
jobName = 'coremantle'
jobDesc = "Hydrostatic equilibrium of a two-layer, fluid planet."
print '\n', jobName.upper(), '-', jobDesc.upper()

# Planet parameters
rPlanet = 1000e3             # Initial guess for outer planet radius (m)
rCore = 500e3                # Initial guess for core radius (m)
matMantle = 'h2oice'         # Mantle material (see <uss>/MATERIALS.md for options)
rhoMantle = 900.             # Initial guess for mantle density (kg/m^3)
matCore = 'h2oice'           # Core material (see <uss>/MATERIALS.md for options)
rhoCore = 2700.              # Initial guess for core density (kg/m^3)
mPlanet = (4.0*pi/3.0) * (rhoCore*rCore**3 + rhoMantle*(rPlanet**3-rCore**3))
rhoPlanet = 3.0*mPlanet/(4.0*pi*rPlanet**3)
gravTime = 1/sqrt(MKS().G*rhoPlanet)
print mPlanet, rhoPlanet, gravTime
