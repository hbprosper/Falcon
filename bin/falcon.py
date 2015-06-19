#!/usr/bin/env python
#-------------------------------------------------------------------------------------
import os, sys, re
from string import *
from math import *
from ROOT import gSystem, TStopwatch, vector
gSystem.Load("libfastjet")
gSystem.Load("libFalcon")
from ROOT import nic
#-------------------------------------------------------------------------------------
PARTONID = {1: "d",
            2: "u",
            3: "s",
            4: "c",
            5: "b",
            11:"e",
            12:"nu_e",
            13:"mu",
            14:"nu_mu",
            15:"tau",
            16:"nu_tau",
            21:"g",
            22:"gamma"}
    
MISSINGETID = {12: "nu_e",
               14: "nu_mu",
               16: "nu_tau"}
def buildMap(filename, events):
    pass

def readEvents(filename, nevents):
    inp = open(filename)
    
    # skip header
    while 1:
        record = split(inp.readline())
        token = record[0]
        if token == '</init>': break
        
    events = []
    nn = 0
    metx = 0.0
    mety = 0.0


    ppx = vector('double')
    ppy = vector('double')
    ppz = vector('double')
    pE  = vector('double')
        
    jetpx = vector('double')
    jetpy = vector('double')
    jetpz = vector('double')
    jetE  = vector('double')

                
    while nn < nevents:
        record = split(inp.readline())
        token = record[0]
        if token == '<event>':
            nn += 1
            if nn % 1000 == 0: print nn
            partons = []
            nonpartons = []
            recobjs = []
            record  = split(inp.readline())
            token = record[0]
            nparticles = atoi(token)
            if nn < 5:
                print "Event: %d" % nn
                print "\tnumber of particles: %d" % nparticles

            # initialize missing ET sums
            metx = 0.0
            mety = 0.0
            continue
        elif token == '</event>':
            
            # missing ET
            energy = sqrt(metx*metx + mety*mety)
            particles = [(12, metx, mety, 0, energy)]
            
            # parton jets
            ppx.clear()
            ppy.clear()
            ppz.clear()
            pE.clear()
            for pid, px, py, pz, E in partons:                
                ppx.push_back(px)
                ppy.push_back(py)
                ppz.push_back(pz)
                pE.push_back(E)
            findJets(vpx, vpy, vpz, vE, jetpx, jetpy, jetpz, jetE)
            for i in xrange(jetpx.size()):
                particles.append( (21, jetpx[i], jetpy[i], jetpz[i], jetE[i]) )

            for p in nonpartons:
                particles.append(p)
            events.append((particles, recobjs))
            continue
        else:
            pid = atoi(token)
            ptype  = record[1]
            px, py, pz, energy, mass = map(atof, record[6:11])
            # truth-level code = 3
            # reco-level code  = 1
            if ptype == '3':
                # exclude beam particles
                if px == py == 0: continue

                # temporary hack: exclude t/W/Z/H
                ID = abs(pid)
                if not PARTONID.has_key(ID): continue
                
                # compute true missing ET
                if MISSINGETID.has_key(ID):
                    metx += px
                    mety += py
                elif (ID < 6)  or (ID == 2):
                    partons.append((pid, px, py, pz, energy))
                else:
                    nonpartons.append((pid, px, py, pz, energy))
            else:
                recobjs.append((pid, px, py, pz, energy))
            if nn < 5:
                if ptype == '3':
                    level = 'truth'
                else:
                    level = 'reco'
                rec = "%6s: %12s %10d\t%10.3f %10.3f %10.3f %10.3f" % \
                  (level, nic.particleName(pid), pid, px, py, pz, energy)
                print rec
    return events
#------------------------------------------------------------------------------------
def main():
    # get mode and LHE file
    argv  = sys.argv[1:]
    build = argv[0][0] in ['b', 'B']
    lhefilename = argv[1]
    if not os.path.exists(lhefilename):
        exit("can't find file %s" % lhefilename)
        
    if len(argv) > 2:
        nevents = atoi(argv[2])
    else:
        nevents = 5000
        
    if build:
        print "=> build map using file %s" % lhefilename
        print "=> map will be written to parton2reco.root"
        events = readEvents(lhefilename, nevents)
        buildMap(lhefilename, events)
    else:
        print "=> simulate events using map in file parton2reco.root"
        print "=> events will be written to %s.root" % lhefilename
        simulateEvents(lhefilename, nevents)
#-------------------------------------------------------------------------------------
try:
    argv = sys.argv[1:]
    if len(argv) < 2:
        exit('''
    Usage:
       falcon.py build|simulate LHE-input-file number-of-events[5000]

       build      map from partons to reco objects
       simulate   use map to simulate reco events
       LHE-input-file
        ''')
    swatch = TStopwatch()
    swatch.Start()
    main()
    print "\telapsed time: %10.2fs" % swatch.RealTime()
    
except KeyboardInterrupt:
    print "ciao!"
    
