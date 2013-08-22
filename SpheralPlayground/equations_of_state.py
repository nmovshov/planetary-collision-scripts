#-------------------------------------------------------------------------------
# Load some of spheral's equation-of-state objects to workspace.
# Run this in iPython for interactive exploration of EOS options in spheral.
#-------------------------------------------------------------------------------
import sys, os
from SolidSpheral3d import *

#-------------------------------------------------------------------------------
# Setup
#-------------------------------------------------------------------------------
# Show signs of life.
print "Loading spheral equations of state..."

# EOS constructors take a units object. I usually work in MKS.
units = PhysicalConstants(1.0,   # Unit length in meters
                          1.0,   # Unit mass in kg
                          1.0)   # Unit time in seconds

#-------------------------------------------------------------------------------
# Tillotson EOS for common materials
#-------------------------------------------------------------------------------
mats = ['Granite', 'Basalt', 'Nylon', 'Pure Ice', '30% Silicate Ice', 'Water']
etamin, etamax = 0.01, 100.0
EOSes = [TillotsonEquationOfState(mat, etamin, etamax, units) for mat in mats]
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
izetl = vector_of_int(1, -1)
initializeANEOS('/proj/nmovshov_hindmost/collisions/ANEOS/ANEOS.INPUT', 'ANEOS.barf', izetl)
SiO2 = ANEOS(0,          # Material number
             1000,       # num rho vals
             1000,       # num T vals
             2000.0,     # minimum density (kg/m^3)
             4000.0,     # maximum density (kg/m^3)
             1.0,        # minimum temperature (K)
             1.0e4,      # maximum temperature (K)
             units)
os.system('rm -f ANEOS.barf')
del izetl

#-------------------------------------------------------------------------------
# A polytropic fluid EOS
#-------------------------------------------------------------------------------
K  = 2e5       # polytropic constant
n  = 1         # polytropic index
mu = 2.2e-3    # mean molecular weight
poly = PolytropicEquationOfStateMKS3d(K,n,mu)
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

