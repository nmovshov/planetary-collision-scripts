#! /proj/nmovshov_hindmost/collisions/SPHERAL/bin/python
#-------------------------------------------------------------------------------
# Compare two methods of producing hydrostatic equilibrium.
#
# This script tests a couple of different ways to initialize node lists inside
# a planet-size body in hydrostatic equilibrium. We don't actually run the SPH
# physics, only generate pressure/density profiles.
#
# To run as an executable check that the shebang line points to the full
# path to spheral's python.
#-------------------------------------------------------------------------------
from math import *
import sys, os
import numpy as np
import mpi # Mike's simplified mpi wrapper
from SolidSpheral3d import *
from GenerateNodeDistribution3d import GenerateNodeDistribution3d
from GenerateNodeDistribution3d import GenerateIcosahedronMatchingProfile3d
from HydroStaticProfile import HydroStaticProfileConstantTemp3d
pcsbase = '' # Edit this with full path to <pcs> if you see an ImportError.
sys.path += ['..',pcsbase,os.getenv('PCSBASE','')]
import shelpers # My module of some helper functions
import PlanetNodeGenerators # My experimental node generators

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

# Diagnostic measures
rho_c = max(planetGenerator.rho)
rho_0 = min(planetGenerator.rho)
P_c = shelpers.pressure(eosPlanet, rho_c, 0)
P_0 = shelpers.pressure(eosPlanet, rho_0, 0)
M_tot = sum(planetGenerator.m)
rvec = np.hypot(planetGenerator.x, np.hypot(planetGenerator.y, planetGenerator.z))
ind = np.argsort(rvec); rvec.sort()
rhovec = np.array(planetGenerator.rho)[ind]
Pvec = np.array([shelpers.pressure(eosPlanet, rho, 0) for rho in rhovec])
r_s = rvec[-1]
rho_s = rhovec[-1]
dPdr = -(Pvec[-1] - Pvec[-2])/(rvec[-1] - rvec[-2])
g = units.G*M_tot/r_s**2

# Print diagnostics
print ""
print "Using quasi-incompressible assumption and inverse EOS we find:"
print "P range   = [{:g} -- {:g}] Pa".format(P_c, P_0)
print "rho range = [{:g} -- {:g}] kg/m^3".format(rho_c, rho_0)
print "surface rho*g - dP/dr = {:g} - {:g} = {:g}".format(
        rho_s*g, dPdr, rho_s*g - dPdr)

# Option two: solve a simplified Lane-Emden style equation, using Cody's class
eostup = (eosPlanet, [0, rPlanet])
stdout = sys.stdout
sys.stdout = open(os.devnull, 'w')
rhoProfile = HydroStaticProfileConstantTemp3d(
                                rho0 = eosPlanet.referenceDensity,
                                rMax = rPlanet,
                                M0 = mPlanet,
                                temp = 200,
                                eostup = eostup,
                                units = units)
sys.stdout = stdout

# Diagnostic measures
rho_c = rhoProfile(0)
rho_0 = rhoProfile(rPlanet)
P_c = shelpers.pressure(eosPlanet, rho_c, 0)
P_0 = shelpers.pressure(eosPlanet, rho_0, 0)
M_tot = mPlanet
rvec = np.array(rhoProfile.soln)[:,0]
rhovec = np.array(rhoProfile.soln)[:,1]
Pvec = np.array([shelpers.pressure(eosPlanet, rho, 0) for rho in rhovec])
r_s = rvec[-1]
rho_s = rhovec[-1]
dPdr = -(Pvec[-1] - Pvec[-2])/(rvec[-1] - rvec[-2])
g = units.G*M_tot/rPlanet**2

# Print diagnostics
print ""
print "Using simplified Lane-Emden integration we find:"
print "P range   = [{:g} -- {:g}] Pa".format(P_c, P_0)
print "rho range = [{:g} -- {:g}] kg/m^3".format(rho_c, rho_0)
print "surface rho*g - dP/dr = {:g} - {:g} = {:g}".format(
        rho_s*g, dPdr, rho_s*g - dPdr)

print ""
