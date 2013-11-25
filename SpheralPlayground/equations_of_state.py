#-------------------------------------------------------------------------------
# Load some of spheral's equation-of-state objects to workspace.
#
# Run this script in python (or preferably iPython) to load all the EOS options
# currently used by uss into the global workspace. There you can interactively
# call the EOS methods and compare different materials. This script can also be
# used to extract the code snippets needed to create EOS objects in spheral runs.
#-------------------------------------------------------------------------------
import sys, os
import SolidSpheral3d as sph # The top-level spheral module importer

#-------------------------------------------------------------------------------
# Setup
#-------------------------------------------------------------------------------
# Show signs of life.
print "Loading spheral equations of state..."

# EOS constructors take a units object. I usually work in MKS.
units = sph.PhysicalConstants(1.0,   # Unit length in meters
                              1.0,   # Unit mass in kg
                              1.0)   # Unit time in seconds

#-------------------------------------------------------------------------------
# Tillotson EOS for common materials
#-------------------------------------------------------------------------------
mats = ['Granite', 'Basalt', 'Nylon', 'Pure Ice', '30% Silicate Ice', 'Water']
etamin, etamax = 0.94, 10.0
pext, pmin, pmax = 0.0, -1e200, 1e200 # these are actually the defaults
EOSes = [sph.TillotsonEquationOfState(mat, 1e-20, 1e20, units,
         etamin_solid = etamin,
         externalPressure = pext, minimumPressure = pmin, maximumPressure = pmax)
         for mat in mats]
granite  = EOSes[0]
basalt   = EOSes[1]
nylon    = EOSes[2]
h2oice   = EOSes[3]
dirtyice = EOSes[4]
water    = EOSes[5]
del EOSes, mats, etamin, etamax

#-------------------------------------------------------------------------------
# M/ANEOS improved SiO2
#-------------------------------------------------------------------------------
izetl = sph.vector_of_int(1, -1)
sph.initializeANEOS('/proj/nmovshov_hindmost/collisions/ANEOS/ANEOS.INPUT',
                                                             'ANEOS.barf', izetl)
etamin, etamax = 0.94, 10
rho0 = 2650
pext, pmin, pmax = 0.0, -1e200, 1e200 # these are actually the defaults
SiO2 = sph.ANEOS(0,           # Material number
                 1000,        # num rho vals
                 1000,        # num T vals
                 etamin*rho0, # minimum density (kg/m^3)
                 etamax*rho0, # maximum density (kg/m^3)
                 1.0,         # minimum temperature (K)
                 1.0e4,       # maximum temperature (K)
                 units,
                 pext, pmin, pmax)
os.system('rm -f ANEOS.barf')
del izetl, etamin, etamax, rho0

#-------------------------------------------------------------------------------
# A polytropic fluid EOS
#-------------------------------------------------------------------------------
K  = 2e5       # polytropic constant
n  = 1         # polytropic index
mu = 2.2e-3    # mean molecular weight
poly = sph.PolytropicEquationOfStateMKS3d(K,n,mu)
del K, n, mu

#-------------------------------------------------------------------------------
# Available materials table
#-------------------------------------------------------------------------------
Materials = {'granite':'Granite solid (Tillotson)',
             'basalt':'Basalt solid (Tillotson)',
             'nylon':'Nylon solid (Tillotson)',
             'h2oice':'Water ice solid (Tillotson)',
             'dirtyice':'30% silicate in water ice (Tillotson)',
             'water':'Liquid water (Tillotson)',
             'SiO2':'Multi phase SiO2 (M/ANEOS)',
             'poly':'Polytrope (n=1 K=2e5)',
            }
print
print "Available materials:"
for k,v in Materials.iteritems():
    print k.ljust(20), v
    pass

