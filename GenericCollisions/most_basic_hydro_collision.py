#! /proj/nmovshov_hindmost/collisions/SPHERAL/bin/python
#-------------------------------------------------------------------------------
# Run two fluid spheres into each other.
# 
# This script launches the most basic type of collision we can imagine. Two
# spherical fluid objects of arbitrary size collide with some specified velocity
# and impact parameter. In spheral terms, the target and impactor are of type
# FluidNodeList, and the only physics package attahced to the integrator is the
# hydro package. Although the target and impactor may select a solid material 
# equation-of-state, they do not possess any elastic strength, and no damage 
# model is attached. No gravity is at play either.
#
# Although this type of collision is not really physically relevant, it may 
# still serve as a sanity check on collisions with more complex physics, and
# also help identify what physical process is most dominant in a given impact.
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
jobName = 'mostbasic'
jobDesc = "Pure hydro collision of fluid, single material spheres."
print '\n', jobName.upper(), '-', jobDesc.upper()

