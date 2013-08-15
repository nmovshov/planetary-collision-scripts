#-------------------------------------------------------------------------------
# Set up interactive python session with spheral modules loaded and basic objects
# initiated.
#
# Run this script in python (or preferably iPython) interpreter to load a minimal
# spheral "session". This helps discover fields and methods of spheral objects,
# interactively, for example by TAB-completion.
# Here is a list of names you can use TAB completion on:
#  eos - the equation of state
#  nodes - the node list
#  generator - the generator
#  db - the node list database 
#  WT - the sph kernel function
#  q - the sph artificial viscosity object
#  hydro - the hydro physics package
#  integrator - the time integrator
#-------------------------------------------------------------------------------
from math import *
import sys, os
import mpi # Mike's simplified mpi wrapper
import shelpers # My module of some helper functions
import SolidSpheral3d as sph # The top-level spheral module importer
from GenerateNodeDistribution3d import GenerateNodeDistribution3d # basic nl-gens
from VoronoiDistributeNodes import distributeNodes3d # the load distributer

#-------------------------------------------------------------------------------
# Construct a minimal spheral simulation structure, consisting of a node list, a
# node lists generator, a node list distributer, a physics package, an integrator,
# and a controller.
#-------------------------------------------------------------------------------
# First, create an equation of state
units = sph.PhysicalConstants(1.0,1.0,1.0)
etamin, etamax = 0.01, 100.0
eos = sph.TillotsonEquationOfState('pure ice',etamin,etamax,units)

# Create an empty node list.
nodes = sph.makeFluidNodeList('nodelist', eos)

# Create a stock generator.
generator = GenerateNodeDistribution3d(20, 20, 20, 
                                       eos.referenceDensity,
                                       distributionType = 'lattice')

# Create a DataBase object to hold the node lists.
db = sph.DataBase()
db.appendNodeList(nodes)

# Create the kernel functions for SPH.
WT = sph.TableKernel(sph.BSplineKernel(), 1000) # one for normal hydro
WTPi = WT                                   # one for artificial viscosity

# Create the artificial viscosity object.
q = sph.MonaghanGingoldViscosity(1.0, 1.0)

# Create the hydro package.
hydro = sph.ASPHHydro(WT,WTPi,q)

# Create the time integrator and attach the physics package to it.
integrator = sph.SynchronousRK2Integrator(db)
integrator.appendPhysicsPackage(hydro)

# Create the controller.
control = sph.SpheralController(integrator, WT)
