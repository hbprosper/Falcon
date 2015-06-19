#ifndef TURTLE_H
#define TURTLE_H
// ---------------------------------------------------------------------------
// Created May 11, 2011 by Harrison Prosper and Sezen Sekmen
// Updated May 21, 2015 HBP - Implement Fill
// ---------------------------------------------------------------------------
#include <vector>
#include <string>
#include <iostream>
#include "TKDTreeBinning.h"
// ---------------------------------------------------------------------------
///
class Turtle
{
public:
  ///
  Turtle();

  ///
  Turtle(std::vector<std::string>& rootfilenames, 
	 std::vector<std::string>& variablenames, 
	 std::string treename,
	 int numberofbins,
	 int numberofpoints=-1);

  virtual ~Turtle();

  /// Histogram data from specified files.
  void Fill(std::vector<std::string>& rootfilenames,
	    std::string weightname="");

  /// Histogram data from specified files.
  void Fill(std::vector<double>& point, double weight=1);
  
  /// Return all bin counts for histogrammed data.
  std::vector<double> GetBinContents() { return _counts; }

  ///
  double GetBinDensity(int bin) { return _btree->GetBinDensity(bin); }

  ///
  double GetBinVolume(int bin)  { return _btree->GetBinVolume(bin); }

  ///
  const double* GetBinCenter(int bin) { return _btree->GetBinCenter(bin); }

  ///
  const double* GetBinWidth(int bin)  { return _btree->GetBinWidth(bin); }

  ///
  void SortBinsByDensity()      { return _btree->SortBinsByDensity(); }
  
  ///
  void Reset();
  
  ClassDef(Turtle,0)
  
 private:
  TKDTreeBinning* _btree;
  std::vector<std::string> _rootfilenames;
  std::vector<std::string> _variablenames; 
  std::string _treename;
  unsigned int _datasize;
  double* _data;
  std::vector<double> _counts;
  int _numberofbins;

  double* _ReadTree(std::vector<std::string>& rootfilenames, 
                    std::vector<std::string>& variablenames, 
                    std::string treename,
                    int numberofbins,
                    int numberofpoints);
};

#endif
