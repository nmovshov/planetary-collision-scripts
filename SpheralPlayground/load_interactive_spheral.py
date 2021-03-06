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
#  control - the spheral controller
#-------------------------------------------------------------------------------
from math import *
import sys, os
import mpi # Mike's simplified mpi wrapper
import SolidSpheral3d as sph # The top-level spheral module importer
from GenerateNodeDistribution3d import GenerateNodeDistribution3d # basic nl-gens
from VoronoiDistributeNodes import distributeNodes3d # the load distributer

pcsbase = '' # Edit this with full path to <pcs> if you see an ImportError.
sys.path += ['..',pcsbase,os.getenv('PCSBASE','')]
import shelpers # My module of some helper functions

#-------------------------------------------------------------------------------
# Construct a minimal spheral simulation structure, consisting of a node list, a
# node lists generator, a node list distributer, a physics package, an integrator,
# and a controller.
#-------------------------------------------------------------------------------
# First, create an equation of state.
units = sph.PhysicalConstants(1.0,1.0,1.0)
eos = shelpers.construct_eos_for_material('h2oice',units)

# Create an empty node list.
nodes = sph.makeFluidNodeList('nodelist', eos)

# Create a stock generator.
generator = GenerateNodeDistribution3d(2, 2, 2, 
                                       eos.referenceDensity,
                                       distributionType = 'lattice')

# Distribute nodes to ranks (suppress with any cl arg to speed things up).
if len(sys.argv) == 1:
    distributeNodes3d((nodes, generator))

# Create a DataBase object to hold the node lists.
db = sph.DataBase()
db.appendNodeList(nodes)

# Create the kernel function for SPH.
WT = sph.TableKernel(sph.BSplineKernel(), 1000)

# Create the artificial viscosity object.
q = sph.MonaghanGingoldViscosity(1.0, 1.0)

# Create the hydro package.
hydro = sph.ASPHHydro(W = WT, Q = q)

# Create the time integrator and attach the physics package to it.
integrator = sph.CheapSynchronousRK2Integrator(db)
integrator.appendPhysicsPackage(hydro)

# Create the controller.
control = sph.SpheralController(integrator, WT)
control.vizBaseName = 'test'
control.vizDir = '.'
