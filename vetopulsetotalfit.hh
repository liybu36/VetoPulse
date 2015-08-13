#ifndef _vetopulsetotalfit_H
#define _vetopulsetotalfit_H

#include <iostream>

#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TChain.h"
#include "TTree.h"
#include "TRint.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"

//#include "vetopulsebasic.hh"
//#include "vetopulsemc.hh"
//#include "vetopulsereal.hh"
#include "vetopulsefitfunc.hh"

using namespace std;

class vetopulsefitfunc;
class vetopulsebasic;

class vetopulsetotalfit:public vetopulsebasic
{
public :
  vetopulsetotalfit():
    MClist(0),Fitlist(0),MCFile(0),MCValue_ntuple(0),
    RealFile(0),FullSpectrum(0),Fit_Total(0),
    startrange(0),endrange(0),primary_npar(0)
  {}

  virtual ~vetopulsetotalfit()
  {}
  
  virtual void Init();
  double Isotope_Sum(int);
  void Isotope_Fraction(int);
  bool Isotope_Activity(int);
  bool Each_Isotope_Sum(int);
  double GetIntegral(double, double); 
  void ScaleSpectrum(TH1F*);
  bool TotalFit(int);
  void SetPrimaryVar();
  void SetFitVar(int);
  double SumFuncs(double*, double*);
  void SaveHistograms();
  bool Load_MCPlots(int);
  bool Load_RealPlots();
  bool BookFitFuncs(int,int);
  void FillFitVars(int);
  double GetLiveTime();
  
  std::vector<vetopulsefitfunc*> func;
  vetopulsefitfunc *total_func;
  
private:
  double startrange, endrange;
  int primary_npar;

  string MCDataFile;
  TFile *MCFile;
  vector<TH1F*> mcplot;
  TFile *RealFile;
  TH1F *FullSpectrum;
  TF1 *Fit_Total;
  TObjArray MClist;
  TObjArray Fitlist;
  TNtuple *MCValue_ntuple;
  vector<TCanvas*> canv;
  vector<TF1*> fitplot;

  double integralsum;
  vector<double> integral_sum; 
  vector<double> isotope_activity;
  vector<double> isotope_fraction;

  vector<double> parvalues;
  vector<string> parnames;
  vector<double> paruplimits;

  ClassDef(vetopulsetotalfit,0);
};

#endif /* _vetopulsetotalfit_ */
