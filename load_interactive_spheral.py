#-------------------------------------------------------------------------------
# Run this in iPython to start an interactive session with spheral modules loaded
#-------------------------------------------------------------------------------
from math import *
import sys, mpi
from SolidSpheral3d import *
from findLastRestart import *
from VoronoiDistributeNodes import distributeNodes3d
from NodeHistory import NodeHistory
from AverageStrain import AverageStrain
