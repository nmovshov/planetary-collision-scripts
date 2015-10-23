import shutil
from SolidSpheral3d import *
from SpheralTestUtilities import *
from findLastRestart import *
from GenerateNodeDistribution3d import *
from math import *
import mpi

import SpheralVoronoiSiloDump
from bisectFunction import bisectFunction

#-------------------------------------------------------------------------------
# Mass function used to create a guassian distribution of nodes
# that are not all equal mass, but where the interacting neighbors
# are very close to equal. This may not be necessary for this problem
# if the companion star does in fact get shredded, but if it doesn't
# having a variable mass distribution saves computing time
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
title("3-D Planet Test")

#-------------------------------------------------------------------------------
# Generic problem params (cgs units)
#-------------------------------------------------------------------------------
commandLine(#distributor = VoronoiDistributeNodes.distributeNodes2d,
            #distributor = DistributeNodes.distributeNodes2d,
            #Geometry
            #rho0    = 2.75,     #g/cc
            rPlanet = 6.0e8,    #cm
            rCore   = 3000.0e5, #cm
            rColl   = 1.0e8,
            vColl   = 3.0e6,    #cm/s
            temp    = 1000.0,
            nrPlanet= 25,
            nrColl  = 20,
            etaMin  = 0.05,
            etaMax  = 5.0,
            #Initial Conditions
            orbiting = False,    # Switch for orbiting/stationary companion

            #Sim params
            useGravity   = True,
            timeStepType = DynamicalTime,
            nPerh        = 2.01,
            CRKSPH         = False,
            gamma        = 5.0/3.0,
            mu           = 1.0,


            #Artificial Viscosity and things
            HydroConstructor = SPHHydro,
            momentumConserving = True,  #for CSPH
            Qconstructor = MonaghanGingoldViscosity,
            Cl = 1.0,
            Cq = 1.0,
            Qhmult = 1.0,
            Qlimiter = False,
            balsaraCorrection = False,
            epsilon2 = 1e-2,
            epsilonTensile = 0.0,
            nTensile = 8,
            negligibleSoundSpeed = 1e-5,
            csMultiplier = 1e-4,
            hmin = 1e-5,
            hmax = 1e20,
            hminratio = 0.05,
            cfl = 0.5,
            XSPH = True,

            #Hydro params
            HEvolution = IdealH,
            densityUpdate = IntegrateDensity,
            #densityUpdate = RigorousSumDensity,
            compatibleEnergyEvolution = True,
            gradhCorrection = False,
            
            plummerLength = 0.1,        # (cm) Plummer softening scale
            opening = 1.0,                 # (dimensionless, OctTreeGravity) opening parameter for tree walk
            fdt = 0.1,                     # (dimensionless, OctTreeGravity) timestep multiplier

            #Time and sim control
            timeStepChoice = AccelerationRatio,
            myIntegrator = CheapSynchronousRK2Integrator3d,
            steps = None,
            goalTime = 200,
            vizTime = 1e-3,
            vizCycle = 100,
            dt = 1e-6,
            dtMin = 1e-9,
            dtMax = 1e5,
            dtGrowth = 2.0,
            dtSample = 1,
            rigorousBoundaries = False,
            verbosedt = False,
            maxSteps = None,
            statsStep = 10,
            redistributeStep = 2000,
            restartStep = 100,
            restoreCycle = None,

            clearDirectories = False,
            dataDir = "earthTest",
            )

#-------------------------------------------------------------------------------
# First convert the command line params to MKS from CGS
#-------------------------------------------------------------------------------
rPlanet = rPlanet * 0.01
rCore   = rCore * 0.01
rColl   = rColl * 0.01
vColl   = vColl * 0.01

#-------------------------------------------------------------------------------
# Construct our base units.
#-------------------------------------------------------------------------------
units = PhysicalConstants(6000000.0,    # Unit length (m)
                          1.0e24,   # Unit mass (kg)
                          1.0)  # Unit time (sec)
#units = CGS()

#-------------------------------------------------------------------------------
# Now convert the cmd params to their unit values
#-------------------------------------------------------------------------------
rPlanet = rPlanet / units.unitLengthMeters
rCore   = rCore / units.unitLengthMeters
rColl   = rColl / units.unitLengthMeters
vColl   = vColl / units.unitLengthMeters
            
dataDir = os.path.join(dataDir,
                       "orbit=%s" % orbiting)


if CRKSPH:
    dataDir = os.path.join(dataDir, "CRKSPH")
    Qconstructor = CRKSPHMonaghanGingoldViscosity
else:
    dataDir = os.path.join(dataDir, "SPH")

restartDir = os.path.join(dataDir, "restarts")
vizDir = os.path.join(dataDir, "visit")
restartBaseName = os.path.join(restartDir, "planetTest")
vizBaseName = "planetTest"


#-------------------------------------------------------------------------------
# Check if the necessary output directories exist.  If not, create them.
#-------------------------------------------------------------------------------
import os, sys
if mpi.rank == 0:
    if clearDirectories and os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    if not os.path.exists(restartDir):
        os.makedirs(restartDir)
    if not os.path.exists(vizDir):
        os.makedirs(vizDir)
mpi.barrier()

#-------------------------------------------------------------------------------
# If we're restarting, find the set of most recent restart files.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    restoreCycle = findLastRestart(restartBaseName)

#-------------------------------------------------------------------------------
# Material properties.
#-------------------------------------------------------------------------------
eosGranite = TillotsonEquationOfState("granite",
                                      etamin=etaMin,
                                      etamax=etaMax,
                                      units=units)
eosIron = TillotsonEquationOfState("iron 130pt",
                                   etamin=etaMin,
                                   etamax=etaMax,
                                   units=units)

rho0 = eosGranite.referenceDensity
rhoC = eosIron.referenceDensity

#-------------------------------------------------------------------------------
# Interpolation kernels.
#-------------------------------------------------------------------------------
WT = TableKernel(BSplineKernel(), 1000)
WTPi = TableKernel(BSplineKernel(), 1000, Qhmult)
output("WT")
output("WTPi")
kernelExtent = WT.kernelExtent

#-------------------------------------------------------------------------------
# Make the NodeList.
#-------------------------------------------------------------------------------
nodesIron = makeFluidNodeList("ironNodes", eosGranite,
                           hmin = hmin,
                           hmax = hmax,
                           hminratio = hminratio,
                           nPerh = nPerh,
                           rhoMin = etaMin*eosIron.referenceDensity,
                           rhoMax = etaMax*eosIron.referenceDensity,
                           topGridCellSize = 1e5)

nodesGranite = makeFluidNodeList("graniteNodes", eosGranite,
                                hmin = hmin,
                                hmax = hmax,
                                hminratio = hminratio,
                                nPerh = nPerh,
                                rhoMin = etaMin*eosGranite.referenceDensity,
                                rhoMax = etaMax*eosGranite.referenceDensity,
                                topGridCellSize = 1e5)

nodeSet = [nodesIron,nodesGranite]
#nodeSet = [nodes1]
for nodes in nodeSet:
    output("nodes.name")
    output("nodes.hmin")
    output("nodes.hmax")
    output("nodes.hminratio")
    output("nodes.nodesPerSmoothingScale")

#-------------------------------------------------------------------------------
# Set the node properties.
#-------------------------------------------------------------------------------
if restoreCycle is None:
    from HydroStaticProfile import HydroStaticProfileConstantTemp3d
    eostup     = (eosIron,[0,rCore],eosGranite,[rCore,rPlanet])
    rhoProfile = HydroStaticProfileConstantTemp3d(rho0,rPlanet,temp,eostup,units)


    for i in xrange(9):
        rr = rPlanet/10.0*i
        print rhoProfile(rr)

    #wait = raw_input("Press Enter to Continue...")
    print "Generating the hydrostatic planet"

    genIron = GenerateIcosahedronMatchingProfile3d(nrPlanet,
                                                     rhoProfile,
                                                     rmin = 0.0,
                                                     rmax = rCore,
                                                     nNodePerh=nPerh,
                                                   rMaxForMassMatching=rPlanet)
    genGranite = GenerateIcosahedronMatchingProfile3d(nrPlanet,
                                                    rhoProfile,
                                                    rmin = rCore,
                                                    rmax = rPlanet,
                                                    nNodePerh=nPerh,
                                                      rMaxForMassMatching=rPlanet)

    msum = mpi.allreduce(sum(genIron.m + [0.0]), mpi.SUM)
    msum += mpi.allreduce(sum(genGranite.m + [0.0]), mpi.SUM)
    assert msum > 0.0
    print "Found planet mass = %g kg." % (msum*units.unitMassKg)

    
    print "Starting node distribution..."
    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes3d
    else:
        from DistributeNodes import distributeNodes3d
    distributeNodes3d((nodesIron,genIron),(nodesGranite,genGranite))
    #distributor((nodes1,genPlanet))
    #distributor((nodes3,genCollider))

    nGlobalNodes = 0
    for n in nodeSet:
        print "Generator info for %s" % n.name
        output("    mpi.allreduce(n.numInternalNodes, mpi.MIN)")
        output("    mpi.allreduce(n.numInternalNodes, mpi.MAX)")
        output("    mpi.allreduce(n.numInternalNodes, mpi.SUM)")
        nGlobalNodes += mpi.allreduce(n.numInternalNodes, mpi.SUM)
    del n
    print "Total number of (internal) nodes in simulation: ", nGlobalNodes

    #wait = raw_input("Press Enter to Continue...")

    # Do some IC stuff here for motion or temperature etc.

    denIron = nodesIron.massDensity()
    denGranite = nodesGranite.massDensity()
    epsIron = nodesIron.specificThermalEnergy()
    epsGranite = nodesGranite.specificThermalEnergy()

    tempIron   = ScalarField("temp", nodesIron)
    tempGranite   = ScalarField("temp", nodesGranite)

    for i in xrange(nodesIron.numInternalNodes):
        tempIron[i] = temp
    for i in xrange(nodesGranite.numInternalNodes):
        tempGranite[i] = temp

    eosIron.setSpecificThermalEnergy(epsIron,denIron,tempIron)
    eosGranite.setSpecificThermalEnergy(epsGranite,denGranite,tempGranite)

'''
    for i in xrange(nodesIron.numInternalNodes):
        epsIron = eosIron.specificThermalEnergy(denIron[i],1000.0)
    for i in xrange(nodesGranite.numInternalNodes):
        epsGranite = eosGranite.specificThermalEnergy(denGranite[i],1000.0)
'''
    
#-------------------------------------------------------------------------------
# Construct a DataBase to hold our node lists.
#-------------------------------------------------------------------------------
db = DataBase()
for n in nodeSet:
    db.appendNodeList(n)
del n
output("db")
output("db.numNodeLists")
output("db.numFluidNodeLists")

#-------------------------------------------------------------------------------
# Construct the artificial viscosities for the problem.
#-------------------------------------------------------------------------------
q = Qconstructor(Cl, Cq)
q.limiter = Qlimiter
q.balsaraShearCorrection = balsaraCorrection
q.epsilon2 = epsilon2
q.negligibleSoundSpeed = negligibleSoundSpeed
q.csMultiplier = csMultiplier
output("q")
output("q.Cl")
output("q.Cq")
output("q.limiter")
output("q.epsilon2")
output("q.negligibleSoundSpeed")
output("q.csMultiplier")
output("q.balsaraShearCorrection")

#-------------------------------------------------------------------------------
# Construct the hydro physics object.
#-------------------------------------------------------------------------------
hydro = HydroConstructor(Q=q,
                         W=WT,
                         WPi=WTPi,
                         cfl = cfl,
                         compatibleEnergyEvolution = compatibleEnergyEvolution,
                         gradhCorrection = gradhCorrection,
                         XSPH = XSPH,
                         densityUpdate = densityUpdate,
                         HUpdate = HEvolution,
                         epsTensile = epsilonTensile,
                         nTensile = nTensile)
output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
output("hydro.gradhCorrection")
output("hydro.XSPH")
output("hydro.densityUpdate")
output("hydro.HEvolution")
output("hydro.epsilonTensile")
output("hydro.nTensile")

packages = [hydro]

#-------------------------------------------------------------------------------
# gravity.
#-------------------------------------------------------------------------------
if useGravity:
    gravity = OctTreeGravity(G = units.G,
                             softeningLength = plummerLength,
                             opening = opening,
                             ftimestep = fdt,
                             timeStepChoice = timeStepChoice)

    packages.append(gravity)

#-------------------------------------------------------------------------------
# Construct a time integrator.
#-------------------------------------------------------------------------------
integrator = myIntegrator(db)
for p in packages:
    integrator.appendPhysicsPackage(p)
integrator.lastDt = dt
if dtMin:
    integrator.dtMin = dtMin
if dtMax:
    integrator.dtMax = dtMax
integrator.dtGrowth = dtGrowth
integrator.rigorousBoundaries = rigorousBoundaries
integrator.verbose = verbosedt
output("integrator")
output("integrator.lastDt")
output("integrator.dtMin")
output("integrator.dtMax")
output("integrator.dtGrowth")
output("integrator.rigorousBoundaries")

#for i in xrange(nodes2.numInternalNodes):
#    print nodes2.Hfield()[i]

#hstats([nodes])
#from SpheralPointmeshSiloDump import dumpPhysicsState
#dumpPhysicsState(integrator, "vizstuff_beforecontroller")

#wait = raw_input("Press Enter to Continue...")
#-------------------------------------------------------------------------------
# Build the controller.
#-------------------------------------------------------------------------------
control = SpheralController(integrator, WT,
                            statsStep = statsStep,
                            restartStep = restartStep,
                            restartBaseName = restartBaseName,
                            vizTime = vizTime,
                            vizDir = vizDir,
                            vizDerivs = True,
                            vizBaseName = "planetTest",
                            restoreCycle = restoreCycle,
                            #restartDir = restartDir,
                            SPH = True)
output("control")

#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
#hstats([nodes])
#dumpPhysicsState(integrator, "vizstuff_aftercontroller")

control.advance(goalTime)

#for i in xrange(nodes2.numInternalNodes):
#    print nodes2.Hfield()[i]

