#ifndef _vetopulsebasic_H
#define _vetopulsebasic_H
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TChain.h"
#include "TTree.h"
#include "TRint.h"
#include "TColor.h"
#include "TNtuple.h"
#include "TString.h"

using namespace std;

class vetopulsebasic:public TObject
{
public :
  vetopulsebasic(){}

  virtual ~vetopulsebasic()
  {}
  
  void   SetMCinputdir(string val) { MCinputdir=val; }
  string GetMCinputdir() { return MCinputdir; }
  void   SetMCFile(string val) { MCFile=val; }
  string GetMCFile() { return MCFile; }
  string GetRealinputdir() { return Realinputdir; }
  void   SetRealinputdir(string val) { Realinputdir=val; }
  void   SetRealFile(string val) { RealFile=val; }
  string GetRealFile() { return RealFile; }
  string Getoutputdir() { return outputdir; }
  void   Setoutputdir(string val) { outputdir=val; }
  void   SetOutFile(string val) { OutFile=val; }
  string GetOutFile() { return OutFile; }

  bool Verifydatafile(TString);
  int Colors(int);
  void SetIsotopes();
  vector<string> GetIsotopes() { return Source; }
  int Search_Isotope(string);
  vector<int> Match_Isotope(vector<string>);
  void Fill_Isotope_Pos(vector<string>);

  vector<double> GetFraction() { return fraction; }
  vector<int> GetNBins() { return NBins; }
  
protected:
  vector<int> Source_Pos;
  vector<string> Source;
  vector<int> PDG;
  vector<int> NBins;
  vector<int> TBins;
  vector<int> Normalization;
  vector<double> fraction;

private:
  string MCFile;
  string MCinputdir;
  string RealFile;
  string Realinputdir;
  string outputdir;
  string OutFile;

  ClassDef(vetopulsebasic,0)
};

inline int vetopulsebasic::Colors(int pos)
{
  vector<int> color;
  color.push_back(TColor::GetColor("#FF2007")); //red         0
  color.push_back(TColor::GetColor("#5A1DE8")); //violet      1
  color.push_back(TColor::GetColor("#000000"));
  color.push_back(TColor::GetColor("#F73CFF")); //pink        5
  color.push_back(TColor::GetColor("#1CFFDF")); //low green   7
  color.push_back(TColor::GetColor("#1485CC")); //blue        4
  color.push_back(TColor::GetColor("#FF791F")); //orange      6
  color.push_back(TColor::GetColor("#AF2FCC")); //dark pink   8
  color.push_back(TColor::GetColor("#E8A60C"));
  color.push_back(TColor::GetColor("#B26618"));
  color.push_back(TColor::GetColor("#79FFFF"));
  color.push_back(TColor::GetColor("#11FF8F"));
  color.push_back(TColor::GetColor("#59FF49")); //green       12    

  return color.at(pos);
}


#endif /* _vetopulsebasic_H */
