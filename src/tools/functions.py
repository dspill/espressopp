# vim:fdm=indent
import espressopp as epp
import mpi4py.MPI as MPI
from espressopp import Int3D, Real3D
import os
import os.path
import re
import time
import glob
import math
import subprocess
from shutil import copy2, rmtree, copyfileobj

def fileOutput(system, integrator, filename, per_atom=True, pressure_tensor=False, full_box=False):
    NPart  = epp.analysis.NPart(system).compute()
    T      = epp.analysis.Temperature(system).compute()
    Pij    = epp.analysis.PressureTensor(system).compute()
    P      = (Pij[0] + Pij[1] + Pij[2]) / 3.0
    step   = integrator.step
    box = system.bc.boxL
    Ek     = (3.0/2.0) * NPart * T
    Etotal = 0.0

    observables = [
        ('step', step),
        ('T', T),
        ('P', P)
    ]
    if pressure_tensor:
        observables += [
            ('P_xx', Pij[0]),
            ('P_yy', Pij[1]),
            ('P_zz', Pij[2]),
            ('P_xy', Pij[3]),
            ('P_xz', Pij[4]),
            ('P_yz', Pij[5])
        ]

    observables.append(('L_x', box[0]))
    if full_box:
            observables.append(('Ly', box[1]))
            observables.append(('Lz', box[2]))

    if per_atom:
        Ek /= NPart
        Etotal += Ek
        observables.append(('E_kin/N', Ek))

        for k in range(system.getNumberOfInteractions()):
            e = system.getInteraction(k).computeEnergy()/NPart
            Etotal += e
            observables.append(('E_pot%d/N' % k, e))

        observables.append(('E_tot/N', Etotal/NPart))

    else:
        Etotal += Ek
        observables.append(('E_kin', Ek))

        for k in range(system.getNumberOfInteractions()):
            e = system.getInteraction(k).computeEnergy()
            Etotal += e
            observables.append(('E_pot%d' % k, e))

        observables.append(('E_tot', Etotal))

    filestream = open(filename, 'a')

    string = ''
    if step == 0:
        string += '# '
        for k in range(len(observables)):
            string += '%-16s ' % ('%d: ' % k + observables[k][0])
    else:
        for k in range(len(observables)):
            string += '%16.9e ' % observables[k][1]

    string += '\n'
    filestream.write(string)


def customWritexyz(filename, system, velocities=True, unfolded=True,
                   append=False):

    if append:
        file = open(filename, 'a')
    else:
        file = open(filename, 'w')

    configurations = epp.analysis.ConfigurationsExt(system)
    configurations.unfolded = unfolded
    configurations.gather()
    configuration = configurations[0]

    if velocities:
        velocities = epp.analysis.Velocities(system)
        velocities.gather()
        velocity = velocities[0]

    numParticles  = int(epp.analysis.NPart(system).compute())
    box_x = system.bc.boxL[0]
    box_y = system.bc.boxL[1]
    box_z = system.bc.boxL[2]
    st = "%d\n%16.9e %16.9e %16.9e\n" % (numParticles, box_x, box_y, box_z)
    file.write(st)

    for pid in configuration:
        xpos   = configuration[pid][0]
        ypos   = configuration[pid][1]
        zpos   = configuration[pid][2]
        if velocities:
            xvel   = velocity[pid][0]
            yvel   = velocity[pid][1]
            zvel   = velocity[pid][2]
            st = "0 %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n" \
                % (xpos, ypos, zpos, xvel, yvel, zvel)
        else:
            st = "0 %16.9e %16.9e %16.9e\n" % (xpos, ypos, zpos)

        file.write(st)


    configurations.clear()
    velocities.clear()
    file.close()


def customWritexyzStream(filestream, system, velocities=True, unfolded=True):
    configurations = epp.analysis.ConfigurationsExt(system)
    configurations.unfolded = unfolded
    configurations.gather()
    configuration = configurations[0]

    if velocities:
        velocities = epp.analysis.Velocities(system)
        velocities.gather()
        velocity = velocities[0]

    numParticles  = int(epp.analysis.NPart(system).compute())
    box_x = system.bc.boxL[0]
    box_y = system.bc.boxL[1]
    box_z = system.bc.boxL[2]
    st = "%d\n%16.9e %16.9e %16.9e\n" % (numParticles, box_x, box_y, box_z)
    filestream.write(st)

    for pid in configuration:
        xpos   = configuration[pid][0]
        ypos   = configuration[pid][1]
        zpos   = configuration[pid][2]
        if velocities:
            xvel   = velocity[pid][0]
            yvel   = velocity[pid][1]
            zvel   = velocity[pid][2]
            st = "0 %16.9e %16.9e %16.9e %16.9e %16.9e %16.9e\n" \
                % (xpos, ypos, zpos, xvel, yvel, zvel)
        else:
            st = "0 %16.9e %16.9e %16.9e\n" % (xpos, ypos, zpos)

        filestream.write(st)

    configurations.clear()
    velocities.clear()
    filestream.flush()


def printInteractions(system):
    for k in range(system.getNumberOfInteractions()):
        print "interaction ", k
        print "    ", system.getInteraction(k)


def arrangeBeadsInSquareShape(system, chain_index, start, monomers_per_chain,
                              bondlen):
    width = math.ceil(math.sqrt(monomers_per_chain))

    x0     = start[0]
    y0     = start[1]
    center = start[2]

    xNew = x0
    yNew = y0
    zNew = center

    dx = bondlen
    dy = 0

    posNew = Real3D(xNew, yNew, zNew)

    for im in range(monomers_per_chain):
        ipart = chain_index * monomers_per_chain + im + 1

        # set new position and type
        system.storage.modifyParticle(ipart, 'pos', posNew)

        dy = 0
        if (im // width) % 2 == 0:
            dx = bondlen
        else:
            dx = -bondlen

        if (im + 1) % width == 0:
            dx = 0
            dy = bondlen

        xNew += dx
        yNew += dy

        posNew = Real3D(xNew, yNew, zNew)


def arrangeChainsInSquareShape(system, num_chains, monomers_per_chain, L, center,
                               bondlen):

    print "Arranging in square shape"
    width = math.ceil(math.sqrt(num_chains))
    separation = L / width

    for ic in range(num_chains):
        i_x = ic % width
        i_y = ic // width

        x0 = i_x * separation
        y0 = i_y * separation

        start = Real3D(x0, y0, center)
        arrangeBeadsInSquareShape(system, ic, start, monomers_per_chain, bondlen)

    system.storage.decompose()


def arrangeBeadsInCubeShape(system, chain_index, start, monomers_per_chain,
                            bondlen):
    width = int(round(math.pow(monomers_per_chain, 1./3.)))

    x0 = start[0]
    y0 = start[1]
    z0 = start[2]

    xNew = x0
    yNew = y0
    zNew = z0

    dx = bondlen
    dy = 0
    dz = 0

    posNew = Real3D(xNew, yNew, zNew)

    for im in range(monomers_per_chain):
        ipart = chain_index * monomers_per_chain + im + 1
        mult = 1

        # set new position and type
        system.storage.modifyParticle(ipart, 'pos', posNew)

        dy = 0
        dz = 0
        # every second row switch directions
        if (im // width) % 2 == 0:
            dx = bondlen
        else:
            dx = -bondlen

        if (im + 1) % width == 0:
            dx = 0
            dy = bondlen

        if (im + 1) % (width * width) == 0:
            mult *= -1
            dx = 0
            dy = 0
            dz = bondlen

        if((im // (width * width)) % 2 == 0):
            mult = 1
        else:
            mult = -1

        xNew += dx
        yNew += mult * dy
        zNew += dz

        posNew = Real3D(xNew, yNew, zNew)


def arrangeChainsInCubeShape(system, num_chains, monomers_per_chain, L,
                             bondlen):

    print "Arranging in cube shape"
    width = math.ceil(math.pow(num_chains, 1./3.))
    separation = float(L) / width

    for ic in range(num_chains):
        i_x = ic % width
        i_y = (ic % (width * width)) // width
        i_z = ic // (width * width)

        x0 = i_x * separation
        y0 = i_y * separation
        z0 = i_z * separation

        start = Real3D(x0, y0, z0)
        arrangeBeadsInCubeShape(system, ic, start, monomers_per_chain, bondlen)

    system.storage.decompose()


def writeParameters(p, filename="info.dat"):
    # write parameters to file
    filestream = open(filename, "a")
    filestream.write(time.strftime("%c") + "\n")
    filestream.write(os.path.basename(__file__) + "\n")

    for key, value in p.iteritems():
        filestream.write("%20s = %s\n" % (key, value))

    filestream.write('\n')
    filestream.close()


def setupSystem(p, xyzfilename=None, phi=0., with_lb=False):
    if xyzfilename is None:
        xyzfilename = "none"

    # Basic setup
    system         = epp.System()
    system.rng     = epp.esutil.RNG()
    system.bc      = epp.bc.OrthorhombicBC(system.rng, p['box'])
    system.skin    = p['skin']

    nodeGrid       = epp.tools.decomp.nodeGrid(epp.MPI.COMM_WORLD.size)
    # nodeGrid       = epp.tools.decomp.nodeGrid(epp.MPI.COMM_WORLD.size,
                                               # p['box'], p['rc'], p['skin'])
    cellGrid       = epp.tools.decomp.cellGrid(p['box'], nodeGrid, p['rc'],
                                               p['skin'])
    system.storage = epp.storage.DomainDecomposition(system, nodeGrid,
                                                     cellGrid)

    # LJcos
    vl          = epp.VerletList(system, p['rc'])
    potential   = epp.interaction.LJcos(phi=phi)
    interaction = epp.interaction.VerletListLJcos(vl)
    interaction.setPotential(0, 0, potential)
    system.addInteraction(interaction)

    integrator     = epp.integrator.VelocityVerlet(system)
    integrator.dt  = p['dt']

    thermostat             = epp.integrator.LangevinThermostat(system)
    thermostat.gamma       = p['langevin_gamma']
    thermostat.temperature = p['temperature']
    integrator.addExtension(thermostat)

    # Configuration
    # Use existing configuration
    if os.path.isfile(xyzfilename):
        print "continuing with configuration " + xyzfilename
        bondlist = readConfiguration(p, system, xyzfilename)
    # Make new configuration
    else:
        print "starting wih NEW configuration"
        bondlist = generateConfiguration(p, system)

    system.storage.decompose()

    # add FENE bonds
    potFENE   = epp.interaction.FENE(K=p['KF'], r0=p['r0'], rMax=p['rMax'])
    interFENE = epp.interaction.FixedPairListFENE(system, bondlist, potFENE)
    system.addInteraction(interFENE)

    # add Lattice Boltzmann
    if with_lb:
        lb = setupLB(p, system, integrator, nodeGrid)
    else:
        lb = None

    return system, integrator, lb


def readConfiguration(p, system, xyzfilename):
    pidf, typef, xposf, yposf, zposf, xvelf, yvelf, zvelf, Lxf, Lyf, Lzf \
        = epp.tools.readxyz(xyzfilename)

    if p['box'] != (Lxf, Lyf, Lzf):
        print("ERROR: input box does not match parameters")
        exit()
        

    if len(pidf) != p['number_of_particles']:
        print("ERROR: input number of particles does not match parameters")
        exit()

    props    = ['id', 'type', 'mass', 'pos', 'v']
    bondlist = epp.FixedPairList(system.storage)
    for i in xrange(p['number_of_chains']):
        chain = []
        bonds = []
        for k in xrange(p['degree_of_polymerization']):
            idx  = i * p['degree_of_polymerization'] + k
            part = [idx, typef[idx], p['mass'],
                    epp.Real3D(xposf[idx], yposf[idx], zposf[idx]),
                    epp.Real3D(xvelf[idx], yvelf[idx], zvelf[idx])]
            chain.append(part)
            if k > 0:
                bonds.append((idx-1, idx))

        system.storage.addParticles(chain, *props)
        system.storage.decompose()
        bondlist.addBonds(bonds)
    return bondlist


def generateConfiguration(p, system):
    props    = ['id', 'type', 'mass', 'pos', 'v']
    bondlist = epp.FixedPairList(system.storage)
    for i in xrange(p['number_of_chains']):
        chain = []
        bonds = []
        for k in xrange(p['degree_of_polymerization']):
            idx  = i * p['degree_of_polymerization'] + k + 1  # TODO +1
            part = [idx, 0, p['mass'], epp.Real3D(0.), epp.Real3D(0.)]
            chain.append(part)

            if k > 0:
                bonds.append((idx-1, idx))

        system.storage.addParticles(chain, *props)
        bondlist.addBonds(bonds)
        system.storage.decompose()

    arrangeChainsInCubeShape(system, p['number_of_chains'],
                             p['degree_of_polymerization'],
                             p['L'], p['bondlen'])
    return bondlist


def setupLB(p, system, integrator, nodeGrid):
    print "adding Lattice-Boltzmann"
    # define a LB grid
    lb = epp.integrator.LatticeBoltzmann(system, Int3D(nodeGrid))
    integrator.addExtension(lb)

    # set parameters of LB fluid (LJ units)
    lb.lbTemp    = p['temperature']  # desired temperature
    lb.nSteps    = p['stepRatio']    # time step contrast (t_{LB} / t_{MD})
    lb.visc_b    = p['visc_b']       # bulk viscosity of LB fluid
    lb.visc_s    = p['visc_s']       # shear viscosity of LB fluid
    lb.profStep  = p['profStep']     # time profiling frequency
    lb.fricCoeff = p['fricCoeff']    # friction coefficient of the coupling

    # initialize populations
    initDen = 1.
    initVel = Real3D(0.)
    initPop = epp.integrator.LBInitPopUniform(system, lb)  # uniform
    initPop.createDenVel(initDen, initVel)

    # screen output
    lboutputScreen = epp.analysis.LBOutputScreen(system, lb)
    outScreen = epp.integrator.ExtAnalyze(lboutputScreen, lb.profStep)
    integrator.addExtension(outScreen)

    # find largest step {{{2
    files = glob.glob('./dump/fluid*.0*')
    if len(files) == 0:
        start_step = 0
    else:
        steps = []
        for file in files:
            new = re.search('.+fluid(.+)\.0\.dat', file)
            steps.append(int(new.group(1)))

        steps.sort()
        start_step = max(steps)

    print "Starting simulation at step " + str(start_step)
    integrator.step = start_step
    return lb


def compactifyDump(filename="./dump.tar", compress=False):
    # find largest step
    files = glob.glob('./dump/fluid*.0*')
    if len(files) == 0:
        start_step = 0
    else:
        steps = []
        for file in files:
            new = re.search('.+fluid(.+)\.0\.dat', file)
            steps.append(int(new.group(1)))

        steps.sort()
        start_step = max(steps)

    print "largest step = " + str(start_step)

    # copy most recent output to separate directory
    if not os.path.isdir("./dump_last"):
        print "creating dump_last"
        os.mkdir("./dump_last")

    files = glob.glob("./dump/*"+str(start_step)+".*")
    print "copying files"
    for file in files:
        copy2(file, "./dump_last/")

    # compress remaining output and remove directory if successful
    if os.path.isdir("./dump"):
        print "compressing dump"
        if compress:
            if subprocess.call(['tar', '-cf', filename, '-I', 'pigz', './dump']) == 0:
                print "compression successful, removing dump"
                # rmtree("./dump")
        else:
            if subprocess.call(['tar', '-cf', filename, './dump']) == 0:
                print "taring successful, removing dump"
                # rmtree("./dump")


def concatenateFiles(outfile, file_list):
    with open(outfile, 'wb') as wfd:
        for f in file_list:
            with open(f, 'rb') as fd:
                copyfileobj(fd, wfd, 1024*1024*10)
                #10MB per writing chunk to avoid reading big file into memory.


def nonInteractingPolymerMelt(num_chains, monomers_per_chain, box=(0, 0, 0),
                              bondlen=0.97, rc=1.5, skin=0.3, dt=0.005, phi=0,
                              KF=30.0, r0=0.0, rMax=1.5, temperature=None,
                              xyzfilename=None, xyzrfilename=None,
                              nodeGrid=None):

    # if given, read configuration from file
    if xyzfilename and xyzrfilename:
        print "ERROR: only one of xyzfilename (only xyz data) or xyzrfilename (additional particle radius data) can be provided."
        sys.exit(1)

    if xyzrfilename:
        pidf, typef, xposf, yposf, zposf, xvelf, yvelf, zvelf, Lxf, Lyf, Lzf, radiusf = epp.tools.readxyzr(xyzrfilename)
        box = (Lxf, Lyf, Lzf)
    elif xyzfilename:
        pidf, typef, xposf, yposf, zposf, xvelf, yvelf, zvelf, Lxf, Lyf, Lzf = epp.tools.readxyz(xyzfilename)
        box = (Lxf, Lyf, Lzf)
    else:
        if box[0] <= 0 or box[1] <= 0 or box[2] <= 0:
            print "WARNING: no valid box size specified, box size set to (100, 100, 100) !"

    system         = epp.System()
    system.rng     = epp.esutil.RNG()
    system.bc      = epp.bc.OrthorhombicBC(system.rng, box)
    system.skin    = skin

    # if no nodeGrid supplied, generate automatically
    if nodeGrid is None:
        nodeGrid       = epp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size)

    cellGrid       = epp.tools.decomp.cellGrid(box, nodeGrid, rc, skin)
    system.storage = epp.storage.DomainDecomposition(system, nodeGrid, cellGrid)

    # add only intrachain interaction!
    interaction = epp.interaction.VerletListLJcos(epp.VerletList(system, cutoff=1.5))
    for k in range(num_chains):
        interaction.setPotential(type1=k, type2=k, potential=epp.interaction.LJcos(phi=phi))
    system.addInteraction(interaction)

    integrator     = epp.integrator.VelocityVerlet(system)
    integrator.dt  = dt
    if(temperature is not None):
        thermostat             = epp.integrator.LangevinThermostat(system)
        thermostat.gamma       = 1.0
        thermostat.temperature = temperature
        integrator.addExtension(thermostat)

    mass     = 1.0

    if xyzfilename:
        props    = ['id', 'type', 'mass', 'pos', 'v']
        bondlist = epp.FixedPairList(system.storage)
        type = 0
        for i in xrange(num_chains):
            chain = []
            bonds = []
            for k in xrange(monomers_per_chain):
                idx  = i * monomers_per_chain + k
                # print idx + 1
                # part = [pidf[idx], type, mass,
                part = [idx + 1, type, mass,
                        epp.Real3D(xposf[idx], yposf[idx], zposf[idx]),
                        epp.Real3D(xvelf[idx], yvelf[idx], zvelf[idx])]
                chain.append(part)
                if k > 0:
                    bonds.append((pidf[idx-1], pidf[idx]))
            type += 1
            system.storage.addParticles(chain, *props)
            system.storage.decompose()
            bondlist.addBonds(bonds)
    else:
        # create new configuration
        props    = ['id', 'type', 'mass', 'pos', 'v']
        vel_zero = epp.Real3D(0.0, 0.0, 0.0)
        bondlist = epp.FixedPairList(system.storage)
        pid      = 1
        type     = 0
        chain    = []
        for i in xrange(num_chains):
            startpos = system.bc.getRandomPos()
            positions, bonds = epp.tools.topology.polymerRW(pid, startpos, monomers_per_chain, bondlen)
            for k in xrange(monomers_per_chain):
                part = [pid + k, type, mass, positions[k], vel_zero]
                chain.append(part)
            pid += monomers_per_chain
            type += 1
            system.storage.addParticles(chain, *props)
            system.storage.decompose()
            chain = []
            bondlist.addBonds(bonds)

    system.storage.decompose()

    # FENE bonds
    potFENE   = epp.interaction.FENE(K=KF, r0=r0, rMax=rMax)
    interFENE = epp.interaction.FixedPairListFENE(system, bondlist, potFENE)
    system.addInteraction(interFENE)

    return system, integrator


def addZConfinement(system, K=30):
    trapPot     = epp.interaction.ZConfinement(K=K)
    interaction = epp.interaction.SingleParticleZConfinement(system, trapPot)
    system.addInteraction(interaction)
