#! /auto/proj/nmovshov_hindmost/collisions/SPHERAL/bin/python
#-------------------------------------------------------------------------------
# A script demostrating the use of spheral's equation-of-state objects.
#-------------------------------------------------------------------------------
from SolidSpheral3d import *
import Gnuplot

#-------------------------------------------------------------------------------
# NAV Show signs of life
#-------------------------------------------------------------------------------
print "Playing around with EOSs..."

#-------------------------------------------------------------------------------
# NAV Build a Tillotson EOS for common materials
#-------------------------------------------------------------------------------
mats = ["Granite", "Basalt", "Nylon", "Pure Ice", "30% Silicate Ice", "Water"]
units = PhysicalConstants(1.0,   # Unit length in meters
                          1.0,   # Unit mass in kg
                          1.0)   # Unit time in seconds
etamin, etamax = 0.01, 100.0
EOSes = [TillotsonEquationOfState(mat, etamin, etamax, units) for mat in mats]

#-------------------------------------------------------------------------------
# NAV Plot pressure-density-temperature triples for granite
#-------------------------------------------------------------------------------
eos = EOSes[0]
rho0 = eos.referenceDensity
rhoMin, rhoMax = 0.9*etamin*rho0, 1.1*etamax*rho0
drho = (rhoMax - rhoMin)/50
rho = [rhoMin + k*drho for k in range(50)]
T = range(200,1001,100)
P = []
for rhok in rho:
	for Tj in T:
		P.append( (rhok, Tj, eos.pressure(rhok,eos.specificThermalEnergy(rhok,Tj))) )

Pplot = Gnuplot.Gnuplot()
Pplot.xlabel("rho (kg/m^3)")
Pplot.ylabel("T (K)")
Pplot(" set grid")
Pdata = Gnuplot.Data(P)
Pplot.splot(Pdata, title="Pressure (Pa)")

#-------------------------------------------------------------------------------
# NAV Plot a P-rho curve for water at reference temperature
#-------------------------------------------------------------------------------
eos = EOSes[5]
T0 = 320
rho0 = eos.referenceDensity
rho = [rho0 + j*rho0/8000 for j in range(101)]
P = []
for rhoj in rho:
	P.append( (rhoj,eos.pressure(rhoj,eos.specificThermalEnergy(rhoj,T0))) )

Prplot = Gnuplot.Gnuplot()
Prplot.xlabel("density (kg/m^3)")
Prplot.ylabel("Pressure (Pa)")
Prplot(" set grid")
Prdata = Gnuplot.Data(P)
Prplot.plot(Prdata)
