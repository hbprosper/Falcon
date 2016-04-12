// ------------------------------------------------------------------------
// FALCON TEST FILE
// ------------------------------------------------------------------------
#include <iostream>
#include <cmath>
#include <cassert>

#include "TKDTreeBinning.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TH2Poly.h"
#include "TStyle.h"
#include "TGraph2D.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TMap.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

#include "FalconTester.h"

using namespace std;

void  FalconTester::Add(string filename)
{
  inputFiles.push_back(filename);
}

// Match Jets with Genjets based on min(dR)
void  FalconTester::DeltaRMatch(double dRcut)
{
  TChain chain("Delphes");
  for(size_t c=0; c < inputFiles.size(); c++)
    {
      cout << "Add: " << inputFiles[c] << endl;
      chain.Add(inputFiles[c].c_str());
    }
  
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  
  for(Int_t entry = 0; entry < numberOfEntries; entry++) //Entry loop
    {
      treeReader->ReadEntry(entry);
      int genjet_entries=branchGenJet->GetEntries();

      double mindR = 999;
      for(int gj=0; gj<genjet_entries; gj++) //GenJet loop
        {
	  Jet* genjet = (Jet*) branchGenJet->At(gj);
	  int jet_entries=branchJet->GetEntries();

	  Jets.push_back(JetMapItem());
	  if ( Jets.size() % 10000 == 0 )
	    cout << "\tjet count = " << Jets.size() << endl;
	  
	  JetMapItem& item = Jets.back(); // NB: get a reference not a copy
	  item.genjet.SetPtEtaPhiM(genjet->PT,
				   genjet->Eta,
				   genjet->Phi,
				   genjet->Mass);
	  
	  for(int j=0; j<jet_entries; j++) //Jet loop
            {
	      Jet *jet = (Jet*) branchJet->At(j);
	      TLorentzVector jet_v;
	      jet_v.SetPtEtaPhiM(jet->PT,jet->Eta,jet->Phi,jet->Mass);
	      double deltar = item.genjet.DeltaR(jet_v);
	      if( deltar < dRcut ) item.jet.push_back(jet_v);
	      if ( deltar < mindR ) mindR = deltar;
            }
	  drmin->Fill(mindR); 
        }
    }
}

void FalconTester::Match(std::string lookfilename,
			 double dRcut)
{
  assert(inputFiles.size()>0);
    
  cout << "1. match genjets to jets..." << endl;
  DeltaRMatch(dRcut);
  
  int totaljets = Jets.size();
  
  std::cout << "\tjet count = " << totaljets << std::endl;

  cout << "2. create lookup table" << endl;
  
  TFile rfile(lookfilename.c_str(), "RECREATE");
  rfile.cd();
  TTree* tree = new TTree("Falcon", "Lookup Table");
  int matched;
  double genPt, genEta, genPhi;
  double Pt, Eta, Phi, Mass;
  tree->Branch("matched", &matched, "matched/I");
  tree->Branch("genPt",   &genPt,   "genPt/D");
  tree->Branch("genEta",  &genEta,  "genEta/D");
  tree->Branch("genPhi",  &genPhi,  "genPhi/D");
  
  tree->Branch("Pt",      &Pt,      "Pt/D");
  tree->Branch("Eta",     &Eta,     "Eta/D");
  tree->Branch("Phi",     &Phi,     "Phi/D");
  tree->Branch("Mass",    &Mass,    "Mass/D");

  for (int entry=0; entry < totaljets; entry++)//loop over entries
    {
      if ( Jets.size() % 10000 == 0 )
	cout << "\tjet count = " << Jets.size() << endl;

      // get references to matched jets
      TLorentzVector&      genjet = Jets[entry].genjet;
      vector<TLorentzVector>& jet = Jets[entry].jet;

      genPt  = genjet.Pt();
      genEta = genjet.Eta();
      genPhi = genjet.Phi();
      
      if ( jet.size() > 0 )
	{
	  matched = 1;
	  Pt  = jet[0].Pt();
	  Eta = jet[0].Eta();
	  Phi = jet[0].Phi();
	  Mass= jet[0].M();
	}
      else
	{
	  matched = 0;
	  Pt  = 0;
	  Eta = 0;
	  Phi = 0;
	  Mass= 0;
	}
	
      rfile.cd();
      tree->Fill();
    }
  cout << "3. write out lookup table" << endl;
  tree->Write();
  rfile.Close();
}


void FalconTester::Build(std::string filename)
{
  cout << "1. open file " << filename << endl;

  TFile rfile(filename.c_str());
  rfile.cd();
  TTree* tree = (TTree*)rfile.Get("Falcon");
  int matched;
  double genPt, genEta, genPhi;
  double Pt, Eta, Phi, Mass;
  tree->SetBranchAddress("matched", &matched);
  tree->SetBranchAddress("genPt",   &genPt);
  tree->SetBranchAddress("genEta",  &genEta);
  tree->SetBranchAddress("genPhi",  &genPhi);
  
  tree->SetBranchAddress("Pt",      &Pt);
  tree->SetBranchAddress("Eta",     &Eta);
  tree->SetBranchAddress("Phi",     &Phi);
  tree->SetBranchAddress("Mass",    &Mass);
  
  int totaljets = tree->GetEntries();
  
  const UInt_t DATASZ  = totaljets;
  const UInt_t DATADIM = 3;
  const UInt_t NBINS   = totaljets;

  cout << "2. allocate ";
  Double_t* smp = new Double_t[DATASZ * DATADIM];
  if(smp==nullptr)
    {
      std::cerr<<"Can not allocate memory: "<<DATASZ * DATADIM<<" doubles\n";
      exit(0);
    }
  cout << sizeof(Double_t)*float(DATASZ*DATADIM)/1000000
       << " Mbytes of memory "
       << endl;

  cout << "3. fill memory with gen level data" << endl;

  for (int entry=0; entry < totaljets; entry++)//loop over entries
    {
      tree->GetEntry(entry);

      smp[DATASZ*0+entry] = genPt;
      smp[DATASZ*1+entry] = genEta;
      smp[DATASZ*2+entry] = genPhi;
    }

  cout << "4. bin gen level data...please be patient!" << endl;

  TStopwatch swatch;
  swatch.Start();
  kdt = new TKDTreeBinning(DATASZ, DATADIM, smp, NBINS);
  cout << "\treal time: " << swatch.RealTime() << " seconds" << endl;
  for (int entry=0; entry < totaljets; entry++)//loop over entries
    {
      tree->GetEntry(entry);
      if ( entry % 20000 == 0 ) cout << "\tjet count = " << entry << endl;
      
      double point[3];
      point[0] = genPt;
      point[1] = genEta;
      point[2] = genPhi;
      int bin  = kdt->FindBin(point);
      RecoJet jet;
      
      if ( matched )
	{
	  jet.PT   = Pt;
	  jet.Eta  = Eta;
	  jet.Phi  = Phi;
	  jet.Mass = Mass;
	}
      else
	{
	  jet.PT   = 0;
	  jet.Eta  = 0;
	  jet.Phi  = 0;
	  jet.Mass = 0;
	}	
      table[bin] = jet; 
    }
  std::cout << "\tjet count = " << totaljets << std::endl;  
  cout << "\tdone!" << endl;
  rfile.Close();
}

RecoJet FalconTester::MapJet(double pt, double eta, double phi)
{
  double point[3] = {pt, eta, phi};
  int index = kdt->FindBin(point);
  if ( table.find(index) != table.end() )
    return table[index];
  else
    return RecoJet();
}
