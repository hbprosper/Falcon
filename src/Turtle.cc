// ---------------------------------------------------------------------------
// File: Turtle.cc
// Description: Bin data in n-dimensions using bins created by recursive
// binary partitioning, such that the count in each bin is the same.
// Created May 11, 2011 by Harrison Prosper and Sezen Sekmen
// Updated May 21, 2015 HBP - Implement Fill
// ---------------------------------------------------------------------------
#include <algorithm>
#include "TMath.h"
#include "TChain.h"
#include "Turtle.h"
// ---------------------------------------------------------------------------
// TKDTreeBinning::TKDTreeBinning()
// {
// }

using namespace std;

Turtle::Turtle()
{
}

Turtle::~Turtle()
{
  if ( _data )  delete _data;
  if ( _btree ) delete _btree;
}

Turtle::Turtle(vector<string>& rootfilenames,
               vector<string>& variablenames,
               string treename,
               int numberofbins,
               int numberofpoints)
  : _btree(0),
    _rootfilenames(rootfilenames),
    _variablenames(variablenames),
    _treename(treename),
    _counts(_variablenames.size(),0)
{
 _btree = new TKDTreeBinning(_datasize,
			     variablenames.size(), 
			     _ReadTree(rootfilenames,
				       variablenames,
				       treename, 
				       numberofbins,
				       numberofpoints), 
			     numberofbins),  
  _numberofbins = _btree->GetNBins();
}

namespace  {
  double zero(double) { return 0; }
};

void Turtle::Fill(std::vector<std::string>& rootfilenames,
		  std::string weightname)
{
  // Histogram data from given files and store values in _counts
  TChain chain(_treename.c_str());
  for (unsigned int i=0; i < rootfilenames.size(); i++)
    chain.Add(rootfilenames[i].c_str());

  int numberofpoints = chain.GetEntries();

  // Allocate enough space for the number of variables
  vector<double> row(_variablenames.size());

  // create a branch for each variable
  for (unsigned int i=0; i < _variablenames.size(); i++)
    chain.SetBranchAddress(_variablenames[i].c_str(), &row[i]);

  double weight = 1.0;
  if ( weightname != "" )
    chain.SetBranchAddress(weightname.c_str(), &weight);
  
  // clear _counts
  transform(_counts.begin(), _counts.end(), _counts.begin(), zero);

  for (int entry=0; entry < numberofpoints; entry++)
    {
      chain.GetEntry(entry);
      if ( entry % 100 == 0 ) 
        cout << "entry: " << entry << endl;
      int bin = _btree->FindBin(&row[0]);
      _counts[bin] += weight;
    }
}

void Turtle::Fill(std::vector<double>& point, double weight)
{
  // Histogram data and store values in _counts
  int bin = _btree->FindBin(&point[0]);
  _counts[bin] += weight;
}

void Turtle::Reset()
{
  // clear _counts
  transform(_counts.begin(), _counts.end(), _counts.begin(), zero);
}

double* Turtle::_ReadTree(vector<string>& rootfilenames, 
                          vector<string>& variablenames, 
                          string treename, 
                          int numberofbins,
                          int numberofpoints)
{
  TChain chain(treename.c_str());
  for (unsigned int i=0; i<rootfilenames.size(); i++)
    chain.Add(rootfilenames[i].c_str());

  _datasize = chain.GetEntries();
  if (numberofpoints > 0)
    _datasize = TMath::Min(_datasize, (unsigned int)numberofpoints); 
  int k = _datasize / numberofbins;
  _datasize = k * numberofbins;

  cout << "datasize: " << _datasize << endl;
  
  // Allocate enough space for the number of points times the 
  // number of variables
  int buffersize = _datasize * variablenames.size();
  _data = new double[buffersize];

  // Allocate enough space for the number of variables
  vector<double> row(variablenames.size());

  for (unsigned int i=0; i < variablenames.size(); i++)
    chain.SetBranchAddress(variablenames[i].c_str(), &row[i]);

  for (unsigned int entry=0; entry < _datasize; entry++)
    {
      chain.GetEntry(entry);
      if ( entry % 100 == 0 ) 
        cout << entry << endl;

      for (unsigned int j=0; j< variablenames.size(); j++)
        _data[entry+j*_datasize] = row[j];
    }
  return &_data[0];
}
