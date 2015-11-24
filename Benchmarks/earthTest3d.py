import shutil
from SolidSpheral3d import *
from SpheralTestUtilities import *
from findLastRestart import *
from GenerateNodeDistribution3d import *
from math import *
import mpi
import harmonic

import SpheralVoronoiSiloDump
from bisectFunction import bisectFunction

#-------------------------------------------------------------------------------
# Mass function used to create a guassian distribution of nodes
# that are not all equal mass, but where the interacting neighbors
# are very close to equal. This may not be necessary for this problem
# if the companion star does in fact get shredded, but if it doesn't
# having a variable mass distribution saves computing time
#-------------------------------------------------------------------------------


class calcHarmonic(object):
    def __init__(self,nodeList,samples,rMax,densityProfile,nominalMass):
        #print "sent in %d %f %f" % (samples,rMax,nominalMass)
        self.nodeList = nodeList
        self.harmonic = harmonic.Harmonic()
        self.samples  = samples;
        self.rMax     = rMax;
        self.time     = []
        self.radii    = []
        self.harmonics= []
    
        self.bins     = []
        self.binwidths= []
        dr = rMax / samples;
        for i in xrange(self.samples):
            r = (i+1.0)*dr
            self.bins.append(r)
            rho = densityProfile(r)
            self.binwidths.append(3.0*pow(nominalMass/rho,1.0/3.0))
            #print "r,dr = %f, %f" %(r,2.0*pow(nominalMass/rho,1.0/3.0))
    def __call__(self,cycle,time,dt):
        self.time.append(time)

        # first get the normalization at each bin
        norm = []
        thisHarm = []
        # (0,0) (1,-1)(1,0)(1,1) (2,-2)(2,-1)(2,0)(2,1)(2,2)
        for i in xrange(self.samples):
            norm.append(0)
            thisHarm.append([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
        for nodes in self.nodeList:
            for i in xrange(nodes.numInternalNodes):
                xi = nodes.positions()[i].x
                yi = nodes.positions()[i].y
                zi = nodes.positions()[i].z
                ri = sqrt(xi*xi + yi*yi + zi*zi)
                vxi= nodes.velocity()[i].x
                vyi= nodes.velocity()[i].y
                vzi= nodes.velocity()[i].z
                vi = (1.0/ri)*(vxi*xi+vyi*yi+vzi*zi)
                for j in xrange(self.samples):
                    rj = self.bins[j]
                    dr = self.binwidths[j]
                    if (abs(ri-rj)<=dr):
                        norm[j] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)
                        thisHarm[j][0] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(0,0,xi,yi,zi)*vi
                        
                        thisHarm[j][1] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(1,-1,xi,yi,zi)*vi
                        thisHarm[j][2] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(1,0,xi,yi,zi)*vi
                        thisHarm[j][3] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(1,1,xi,yi,zi)*vi
                        
                        thisHarm[j][4] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(2,-2,xi,yi,zi)*vi
                        thisHarm[j][5] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(2,-1,xi,yi,zi)*vi
                        thisHarm[j][6] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(2,0,xi,yi,zi)*vi
                        thisHarm[j][7] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(2,1,xi,yi,zi)*vi
                        thisHarm[j][8] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(2,2,xi,yi,zi)*vi
                        
                        thisHarm[j][9] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(3,-3,xi,yi,zi)*vi
                        thisHarm[j][10] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(3,-2,xi,yi,zi)*vi
                        thisHarm[j][11] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(3,-1,xi,yi,zi)*vi
                        thisHarm[j][12] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(3,0,xi,yi,zi)*vi
                        thisHarm[j][13] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(3,1,xi,yi,zi)*vi
                        thisHarm[j][14] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(3,2,xi,yi,zi)*vi
                        thisHarm[j][15] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(3,3,xi,yi,zi)*vi
                        
                        thisHarm[j][16] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(4,-4,xi,yi,zi)*vi
                        thisHarm[j][17] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(4,-3,xi,yi,zi)*vi
                        thisHarm[j][18] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(4,-2,xi,yi,zi)*vi
                        thisHarm[j][19] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(4,-1,xi,yi,zi)*vi
                        thisHarm[j][20] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(4,0,xi,yi,zi)*vi
                        thisHarm[j][21] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(4,1,xi,yi,zi)*vi
                        thisHarm[j][22] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(4,2,xi,yi,zi)*vi
                        thisHarm[j][23] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(4,3,xi,yi,zi)*vi
                        thisHarm[j][24] += (1.0/ri)*pow(nodes.mass()[i]/nodes.massDensity()[i],2.0/3.0)*self.harmonic.Y(4,4,xi,yi,zi)*vi
                        break
            
        for i in xrange(self.samples):
            if(norm[i]>0):
                norm[i] = 4.0*pi/norm[i]
            for j in xrange(25):
                thisHarm[i][j] = abs(thisHarm[i][j])*norm[i]
        
        collapse = []
        for j in xrange(self.samples):
            collapse.append([thisHarm[j][0],
                             thisHarm[j][1]+thisHarm[j][2]+thisHarm[j][3],
                             thisHarm[j][4]+thisHarm[j][5]+thisHarm[j][6]+thisHarm[j][7]+thisHarm[j][8],
                             thisHarm[j][9]+thisHarm[j][10]+thisHarm[j][11]+thisHarm[j][12]+thisHarm[j][13]+thisHarm[j][14]+thisHarm[j][15],
                             thisHarm[j][16]+thisHarm[j][17]+thisHarm[j][18]+thisHarm[j][19]+thisHarm[j][20]+thisHarm[j][21]+thisHarm[j][22]+thisHarm[j][23]+thisHarm[j][24]])
            
        print self.bins
        print collapse
        return


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
            mPlanet = 5.97e27,  # g
            rColl   = 1.0e8,
            vColl   = 3.0e6,    #cm/s
            temp    = 1000.0,
            nrPlanet= 30,
            nrColl  = 20,
            etaMin  = 0.05,
            etaMax  = 50.0,
            #Initial Conditions
            orbiting = False,    # Switch for orbiting/stationary companion

            #Sim params
            useGravity   = True,
            timeStepType = DynamicalTime,
            nPerh        = 1.51,
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
            dataDir = "earthTest3d",
            historyFileName = "earthHistory.txt",
            #stuff for harmonic oscillation test
            sampleBins = 5,
            checkEvery = 100,
            harmonicFile = "harmonics.dat",
            
            )

#-------------------------------------------------------------------------------
# First convert the command line params to MKS from CGS
#-------------------------------------------------------------------------------
rPlanet = rPlanet * 0.01
mPlanet = mPlanet * 0.001
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
mPlanet = mPlanet / units.unitMassKg
rCore   = rCore / units.unitLengthMeters
rColl   = rColl / units.unitLengthMeters
vColl   = vColl / units.unitLengthMeters
            
dataDir = os.path.join(dataDir,"n=%d" % nrPlanet)
dataDir = os.path.join(dataDir,"R=%3.2f,M=%3.2f" % (rPlanet,mPlanet))


if CRKSPH:
    dataDir = os.path.join(dataDir, "CRKSPH")
    Qconstructor = CRKSPHMonaghanGingoldViscosity
    HydroConstructor = CRKSPHHydro
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
eosGranite = TillotsonEquationOfState("basalt",
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
WT = TableKernel(NBSplineKernel(5), 1000)
WTPi = TableKernel(NBSplineKernel(5), 1000, Qhmult)
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
    from HydroStaticProfile import EarthLikeProfileConstantTemp3d
    eostup      = (eosIron,[0,rCore],eosGranite,[rCore,rPlanet])
    rhoProfile  = EarthLikeProfileConstantTemp3d(rho0,rPlanet,mPlanet,temp,eostup,units)
    print "Found new rMax = {0:3.3e} Old was {1:3.3e}".format(rhoProfile.rMax,rPlanet)
    rPlanet     = rhoProfile.rMax


    for i in xrange(100):
        rr = rPlanet/100.0*i
        print "%e %e" %(rr,rhoProfile(rr)*(units.unitMassKg/(units.unitLengthMeters)**3*0.001))

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

    massIron = nodesIron.mass()
    massGrainte = nodesGranite.mass()
    denIron = nodesIron.massDensity()
    denGranite = nodesGranite.massDensity()
    epsIron = nodesIron.specificThermalEnergy()
    epsGranite = nodesGranite.specificThermalEnergy()

    tempIron   = ScalarField("temp", nodesIron)
    tempGranite   = ScalarField("temp", nodesGranite)

    mScale = mPlanet/msum

    for i in xrange(nodesIron.numInternalNodes):
        massIron[i] = massIron[i]*mScale
        tempIron[i] = temp
    for i in xrange(nodesGranite.numInternalNodes):
        massGrainte[i] = massGrainte[i]*mScale
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
                         #gradhCorrection = gradhCorrection,
                         XSPH = XSPH,
                         densityUpdate = densityUpdate,
                         HUpdate = HEvolution,
                         #epsTensile = epsilonTensile,
                         #nTensile = nTensile
                         )
output("hydro")
output("hydro.kernel()")
output("hydro.PiKernel()")
output("hydro.cfl")
output("hydro.compatibleEnergyEvolution")
#output("hydro.gradhCorrection")
output("hydro.XSPH")
output("hydro.densityUpdate")
output("hydro.HEvolution")
#output("hydro.epsilonTensile")
#output("hydro.nTensile")

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
                            vizStep=vizCycle,
                            vizDir = vizDir,
                            vizDerivs = True,
                            vizBaseName = "planetTest",
                            restoreCycle = restoreCycle,
                            #restartDir = restartDir,
                            SPH = True)
output("control")

harmonicCalculator = calcHarmonic(nodeSet,sampleBins,rPlanet,rhoProfile,genGranite.m0)

if (checkEvery is not None):
    control.appendPeriodicTimeWork(harmonicCalculator,checkEvery)


#-------------------------------------------------------------------------------
# Advance to the end time.
#-------------------------------------------------------------------------------
#hstats([nodes])
#dumpPhysicsState(integrator, "vizstuff_aftercontroller")

control.advance(goalTime)

control.step(1)
control.doPeriodicWork(force=True)


control.advance(goalTime)
control.doPeriodicWork(force=True)

rank = mpi.rank
if (checkEvery is not None and rank ==0):
    f = open(dataDir + "/" + harmonicFile,'w')
    for i in xrange(len(harmonicCalculator.time)):
        f.write("{0}\n".format(harmonicCalculator.time[i]))
    f.close()

control.conserve.writeHistory(historyFileName)

#for i in xrange(nodes2.numInternalNodes):
#    print nodes2.Hfield()[i]

