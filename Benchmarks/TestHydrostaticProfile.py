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
nxPlanet = 40                # Nodes across a diameter
mPlanet = rhoPlanet*4.0*pi/3.0*rPlanet**3
gravTime = 1/sqrt(MKS().G*rhoPlanet)

# Equation of state
units = PhysicalConstants(1.0, # unit length in meters
                          1.0, # unit mass in kilograms
                          1.0) # unit time in seconds
eosPlanet = shelpers.construct_eos_for_material(matPlanet,units)
assert eosPlanet.valid()
# Optionally, provide non-default values to the following
#eosPlanet.etamin_solid = 0.94 # default is 0.94
#eosPlanet.minimumPressure = 0.0 # default is -1e+200

# Option one: calculate a pressure/density profile using my hydrostaticizing
# functions on an HCP generator.
planetGenerator = PlanetNodeGenerators.HexagonalClosePacking(
                    nx = nxPlanet,
                    rho = rhoPlanet,
                    scale = 2*rPlanet,
                    rMin = 0.0,
                    rMax = rPlanet)
planetGenerator.EOS = eosPlanet
shelpers.hydrostaticize_one_layer_planet(planetGenerator)
rho_c = max(planetGenerator.rho)
P_c = shelpers.pressure(eosPlanet, rho_c, 0)

print ""
print "Using quasi-incompressible assumption and inverse EOS we find:"
print "P_center = {:g} Pa; rho_center = {:g} kg/m^3".format(P_c, rho_c)
