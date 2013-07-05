#! /auto/proj/nmovshov_hindmost/collisions/SPHERAL/bin/python
#-------------------------------------------------------------------------------
# A script demostrating the use of spheral's equation-of-state objects.
# Run this in iPython for interactive exploration of EOS options in spheral.
#-------------------------------------------------------------------------------
from SolidSpheral3d import *
import Gnuplot

#-------------------------------------------------------------------------------
# NAV Show signs of life
#-------------------------------------------------------------------------------
print "Loading spheral equations of state..."

#-------------------------------------------------------------------------------
# NAV Build a Tillotson EOS for common materials
#-------------------------------------------------------------------------------
mats = ["Granite", "Basalt", "Nylon", "Pure Ice", "30% Silicate Ice", "Water"]
units = PhysicalConstants(1.0,   # Unit length in meters
                          1.0,   # Unit mass in kg
                          1.0)   # Unit time in seconds
etamin, etamax = 0.01, 100.0
EOSes = [TillotsonEquationOfState(mat, etamin, etamax, units) for mat in mats]
granite  = EOSes[0]
basalt   = EOSes[1]
nylon    = EOSes[2]
h2oice   = EOSes[3]
dirtyice = EOSes[4]
water    = EOSes[5]
del EOSes

#-------------------------------------------------------------------------------
# NAV Print available materials table
#-------------------------------------------------------------------------------
Materials = {'granite':'Granite solid (Tillotson)',
             'basalt':'Basalt solid (Tillotson)',
	     'nylon':'Nylon solid (Tillotson)',
	     'h2oice':'Water ice solid (Tillotson)',
	     'dirtyice':'30% silicate in water ice (Tillotson)',
	     'water':'Liquid water (Tillotson',
	     }
print
print "Available materials:"
for k,v in Materials.iteritems():
	print k.ljust(20), v


