#ifndef _vetopulsefitfunc_H
#define _vetopulsefitfunc_H
#include <map>

#include "TH1F.h"
#include "TH2F.h"
#include "TMinuit.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TChain.h"
#include "TTree.h"
#include "TRint.h"
#include "TColor.h"
#include "TNtuple.h"

#include "vetopulsebasic.hh"

using namespace std;

class vetopulsebasic;

class vetopulsefitfunc:public vetopulsebasic
{
public :
  vetopulsefitfunc(){}

  vetopulsefitfunc(TH1F* f, int k): fMCFunc(f), npar(k)
  {}
  
  virtual ~vetopulsefitfunc()
  {}
  
  double scint_e_quenching(double, double*);
  double response_function(double, double, double*,int);
  
  double Eval(double*,double*);
  double FitC14(double*,double*);
  double FitFunc(double*,double*);
  double FitFunc_1st(double*,double*);
  double FitFunc_2nd(double*,double*);
  //  double SumFuncs(double*,double*);
  //  void SetFitPlot(vector<TF1*> var) { fitplot=var; }
  
private:
  TH1F* fMCFunc;
  int   npar;  
  //  vector<TF1*> fitplot;  

  ClassDef(vetopulsefitfunc,0);
};

#endif /* _vetopulsefitfunc_ */


  /*
  double Co60Fit(double*,double*);
  double Co57Fit(double*,double*);
  double K40Fit(double*,double*);
  double Th232Fit(double*,double*);
  double Th232LowerFit(double*,double*);
  double U235Fit(double*,double*);
  double U235LowerFit(double*,double*);
  double U238Fit(double*,double*);
  double U238LowerFit(double*,double*);
  double TotalFit(double*,double*);
  */  

