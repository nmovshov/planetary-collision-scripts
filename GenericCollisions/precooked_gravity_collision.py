#! /proj/nmovshov_hindmost/collisions/SPHERAL/bin/python
#-------------------------------------------------------------------------------
# first stab at loading pre-built planets. this is a rushed job to be fixed after
# agu abstract is due
#-------------------------------------------------------------------------------
from math import *
import sys, os
import random
import mpi # Mike's simplified mpi wrapper
import shelpers # My module of some helper functions
from SolidSpheral3d import *
from findLastRestart import *
from VoronoiDistributeNodes import distributeNodes3d
from NodeHistory import NodeHistory
from GenerateNodeDistribution3d import GenerateNodeDistribution3d

jobName = 'precooked'
jobDesc = "load and run precooked planet(s)."
print '\n', jobName.upper(), '-', jobDesc.upper()


