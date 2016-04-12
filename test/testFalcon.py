#!/usr/bin/env python
# ----------------------------------------------------------------------------
import os, sys
from sys import exit
from histutil import setStyle, mkhist1
from time import sleep
from ROOT import *
# ----------------------------------------------------------------------------
print gSystem.Load("libFalcon")
# ----------------------------------------------------------------------------
def build(falcon, inputfiles):
    for filename in inputfiles:
        if not os.path.exists(filename):
            exit("can't open input file %s" % filename)        
        falcon.Add(filename)
    # match partons to reco objects and
    # save mapping.
    falcon.Match()
# ----------------------------------------------------------------------------
def simulate(falcon, inputfiles):
    code = '''#include <string>
#include <cassert>
#include "TChain.h"    
#include "TClonesArray.h"    
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
using namespace std;

struct Reader : public ExRootTreeReader
{
  Reader() : chain(new TChain("Delphes")),
  reader(0), b_genjet(0), b_jet(0), numberOfEntries(0)
  {
    reader   = new ExRootTreeReader(chain);
    assert(reader!=0);
  }
  ~Reader() { delete reader; delete chain; }

  void Add(string filename)
  {
    chain->Add(filename.c_str());
    if ( !b_jet )
     { 
       b_jet    = reader->UseBranch("Jet");
       assert(b_jet!=0);
     }

    if ( !b_genjet )
     {
       b_genjet = reader->UseBranch("GenJet");
       assert(b_genjet!=0);
     }

    if ( numberOfEntries <= 0 )
      numberOfEntries = reader->GetEntries();        
  }

  void Read(int entry)
  {
    reader->ReadEntry(entry);
    numberOfGenJets = b_genjet->GetEntries();
    numberOfJets    = b_jet->GetEntries();
  }

  Jet* GetGenJet(int entry)
  {
    return (Jet*)b_genjet->At(entry);
  }

  Jet* GetJet(int entry)
  {
    return (Jet*)b_genjet->At(entry);
  }


  TChain* chain;
  ExRootTreeReader* reader;
  TClonesArray* b_genjet;
  TClonesArray* b_jet;
  
  int numberOfEntries;
  int numberOfGenJets;
  int numberOfJets;
};
    '''
    gROOT.ProcessLine(code)

    reader = Reader()
    for filename in inputfiles:
        print "=> Add: %s" % filename
        reader.Add(filename)
    print "\n=> numberOfEntries: %d\n" % reader.numberOfEntries
     
    # for now, we have to build the k-d tree
    # every time we want to use it. But, once we
    # can store k-d trees in Root files this won't
    # necessary and the Build() will be done in Match()
    falcon.Build()

    # ----------------------------------------------------------------------
    
    setStyle()
    
    rfile = TFile("histograms.root", "RECREATE")
    
    # transverse momenta
    fjet1pt  = mkhist1('fjet1pt',  '#font[12]{p}_{T,jet1}','', 50, 0, 2000)
    fjet2pt  = mkhist1('fjet2pt',  '#font[12]{p}_{T,jet2}','', 50, 0, 2000)
    fjet3pt  = mkhist1('fjet3pt',  '#font[12]{p}_{T,jet3}','', 50, 0, 1000)

    hjet1pt  = mkhist1('hjet1pt',  '#font[12]{p}_{T,jet1}','', 50, 0, 2000)
    hjet2pt  = mkhist1('hjet2pt',  '#font[12]{p}_{T,jet2}','', 50, 0, 2000)
    hjet3pt  = mkhist1('hjet3pt',  '#font[12]{p}_{T,jet3}','', 50, 0, 1000)


    for entry in xrange(reader.numberOfEntries):
        reader.Read(entry)
        if entry % 1000 == 0: print entry

        # Falcon
        fjet = []
        for ii in xrange(reader.numberOfGenJets):
            genjet = reader.GetGenJet(ii)            
            jet = falcon.MapJet(genjet.PT, genjet.Eta, genjet.Phi)
            if not ( jet.PT > 30 ): continue
            if not ( abs(jet.Eta) < 4.5 ): continue

            fjet.append(jet)

        if len(fjet) > 2:
            fjet1pt.Fill(fjet[0].PT)
            fjet2pt.Fill(fjet[1].PT)
            fjet3pt.Fill(fjet[2].PT)    
     

        # Delphes 
        rjet = []
        for ii in xrange(reader.numberOfJets):
            jet = reader.GetJet(ii)
            if jet == None: continue # THIS SHOULD NOT HAPPEN!!!

            if not ( jet.PT > 30 ): continue
            if not ( abs(jet.Eta) < 4.5 ): continue

            rjet.append(jet)

        if len(rjet) > 2:
            hjet1pt.Fill(rjet[0].PT)
            hjet2pt.Fill(rjet[1].PT)
            hjet3pt.Fill(rjet[2].PT)    
     
    # plot
    
    hist = []
    for ii, ymax in enumerate([250, 350, 700]):
        jj = ii + 1
        
        h = eval('hjet%dpt' % jj)
        h.SetLineWidth(2)
        h.SetMarkerSize(0.8)
        h.SetMaximum(ymax)

        f = eval('fjet%dpt' % jj)                
        f.SetLineWidth(3)
        f.SetLineColor(kBlue)

        name = 'falcon_fig%d' % jj
        c = TCanvas(name, name, 10+500*ii, 10+50*ii, 500, 500)    
        c.cd()
        h.Draw('ep')
        f.Draw('hist same')
        c.Update()
        c.SaveAs('.pdf')
        hist.append((c, h, f))
    sleep(10)
            
    rfile.Close()
    
# ----------------------------------------------------------------------------
def main():
    falcon = FalconTester()
    
    if len(sys.argv) > 1:
        simulate(falcon,
                 ['../data/H2.root'])        
    else:        
        build(falcon,
              ["../data/ttbar13TeV.root",
               "../data/H213TeV2.root"])

# ---------------------------------------------------------------------------
try:
    main()
except KeyboardInterrupt:
    print "bye!"

