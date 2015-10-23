#! /proj/nmovshov_hindmost/collisions/SPHERAL/bin/python
#-------------------------------------------------------------------------------
# Compare different methods of producing hydrostatic equilibrium.
#
# This script tests a couple of different ways to initialize node lists inside
# a planet-size body in hydrostatic equilibrium. We don't actually run the SPH
# physics, only generate initial node lists.
#
# To run as an executable check that the shebang line points to the full
# path to spheral's python.
#-------------------------------------------------------------------------------
from math import *
import sys, os
import mpi # Mike's simplified mpi wrapper
from SolidSpheral3d import *
from GenerateNodeDistribution3d import GenerateNodeDistribution3d
from GenerateNodeDistribution3d import GenerateIcosahedronMatchingProfile3d
pcsbase = '' # Edit this with full path to <pcs> if you see an ImportError.
sys.path += ['..',pcsbase,os.getenv('PCSBASE','')]
import shelpers # My module of some helper functions
import PlanetNodeGenerators # New experimental node generators

# Say something
print "TEST OF HYDROSTATIC DENSITY PROFILES"

# Planet parameters
rPlanet = 1000e3             # Initial guess for planet radius (m)
rhoPlanet = 2700.            # Initial guess for planet density (kg/m^3)
matPlanet = 'basalt'         # Planet material (see <pcs>/MATERIALS.md for options)
mPlanet = rhoPlanet*4.0*pi/3.0*rPlanet**3
gravTime = 1/sqrt(MKS().G*rhoPlanet)
