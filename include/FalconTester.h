#ifndef FALCON_TESTER_H
#define FALCON_TESTER_H
// ---------------------------------------------------------------------------
// Created March 2016 by Sergei Gleyzer
// ---------------------------------------------------------------------------
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include "TObject.h"
#include "TH1D.h"
#include "TLorentzVector.h"
#include "TKDTreeBinning.h"
#include "classes/DelphesClasses.h"
// ---------------------------------------------------------------------------

struct RecoJet
{
RecoJet() : PT(0), Eta(0), Phi(0), Mass(0) {}
  double PT;
  double Eta;
  double Phi;
  double Mass;
};

/// Map (GenJet, Jet, mindr).
struct JetMapItem
{
  JetMapItem()
  : genjet(TLorentzVector()), jet(std::vector<TLorentzVector>())
  {}
  ~JetMapItem() {}
  
  TLorentzVector genjet;
  std::vector<TLorentzVector> jet;
};

struct FalconTester
{
  FalconTester()
  : drmin(new TH1D("drmin", "#DeltaR", 100, 0, 1)),
    Jets(std::vector<JetMapItem>()),
    kdt(0),
    table(std::map<int, RecoJet>()),
    inputFiles(std::vector<std::string>()){}
  ~FalconTester() {}

  void  Add(std::string filename);

  void  Match(std::string lookfilename="JetLookupTable.root",
	      double dRcut=0.35);

  void  Build(std::string lookfilename="JetLookupTable.root");

  void  DeltaRMatch(double dRcut);

  RecoJet MapJet(double  pt, double  eta, double  phi);
  
  TH1D* drmin; 
  std::vector<JetMapItem> Jets;
  TKDTreeBinning* kdt;
  std::map<int, RecoJet> table;
  std::vector<std::string> inputFiles;
};

#endif
