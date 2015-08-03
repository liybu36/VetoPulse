#include "vetopulsefitfunc.hh"

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

#include "TH1F.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TRint.h"
#include "TMath.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TPad.h"
#include "TNtuple.h"
#include "TString.h"

using namespace std;

double vetopulsefitfunc::scint_e_quenching (double e_keV, double* params) {
  double birks_value[7][6] = {{0.49041,0.21600,0.10725,-1.46804e-3,0.22468,0.09544}, //0.006
			      {0.42149,0.18731,0.09624,-0.52346e-3,0.19977,0.08311}, //0.008
			      {0.36851,0.16398,0.08816,-0.45562e-3,0.17582,0.07431}, //0.010
			      {0.32903,0.14404,0.08059,-0.04536e-3,0.15623,0.06611}, //0.012
			      {0.29668,0.12872,0.07477,0.15992e-3,0.14091,0.05978}, //0.014
			      {0.27020,0.11663,0.06985,0.34938e-3,0.12933,0.05453}, //0.016
			      {0.24808,0.10646,0.06576,0.46844e-3,0.11938,0.04998}  //0.018
  };
  int k=0;
  double A1 = birks_value[k][0];
  double A2 = birks_value[k][1];
  double A3 = birks_value[k][2];
  double A4 = birks_value[k][3];
  double A5 = birks_value[k][4];
  double A6 = birks_value[k][5];
  /*  
  double kB =0;// params[1];
  double A1 = -0.6292-0.2181*log(kB);
  double A2 = -0.3057-0.1024*log(kB);
  double A3 = -0.0673-0.03353*log(kB);
  double A4 = 0.009876+0.002276*log(kB);
  double A5 = -0.2814-0.09964*log(kB);
  double A6 = -0.09962-0.0376*log(kB);
  */
  return ( (A1 + A2*log(e_keV) + A3*log(e_keV)*log(e_keV) +
	    A4*log(e_keV)*log(e_keV)*log(e_keV))/(1 + A5*log(e_keV)
						  + A6*log(e_keV)*log(e_keV) + A4*log(e_keV)*log(e_keV)*log(e_keV)) );
}

double vetopulsefitfunc::response_function (double q, double energy, double* params, int choice) {
  double ly_mean = params[1];
  double spe_var = params[2];
  double ly_var = params[3];
  double baseline_mean = params[5];
  double baseline_var = params[6];
  double threshold = params[7];
  double pe;
  double gaus_mean;

  pe = energy*ly_mean;
  //  if (energy == 0) pe = energy*ly_mean;
  //pe = energy*ly_mean*scint_e_quenching(energy, params);

  if (choice == 1) gaus_mean = pe + (baseline_mean + threshold);
  else gaus_mean = pe + baseline_mean;

  double gaus_var = baseline_var + (1 + spe_var)*pe + ly_var*pe*pe;
  //  double gaus_var = baseline_var*(1+spe_var)*gaus_mean+ly_var*TMath::Power(gaus_mean,2);
  double gaus_var_inv = 1.0 / gaus_var;
  double arg = -0.5*(q - gaus_mean)*(q - gaus_mean)*gaus_var_inv;

  return 0.3989422804 * sqrt(gaus_var_inv) * exp(arg);
}

double vetopulsefitfunc::Eval(double *x, double *params)
{
  int nBins = fMCFunc->GetNbinsX();
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy   = fMCFunc->GetBinCenter(i);
      spectrum = fMCFunc->GetBinContent(i);
      result  += spectrum*response_function(q,energy,params,choice);
    }
  double rate = params[npar];
  return result;
}

double vetopulsefitfunc::FitFunc(double *x, double *params)
{
  double result = Eval(x, params);
  double delta = fMCFunc->GetBinWidth(1);
  double rate   = params[npar];
  return result*(rate)*params[0]*delta;    
}

double vetopulsefitfunc::FitC14(double *x,double *params)
{
  int nBins = fMCFunc->GetNbinsX();
  double delta = fMCFunc->GetBinWidth(1);
  int choice = 1;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = fMCFunc->GetBinCenter(i);
      spectrum = fMCFunc->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double C14Rate = params[npar];

  return result*(C14Rate*(2100.e-9)/449177)*params[0]*delta + params[4]; //fieldon
  //  return result*(C14Rate*(2100.e-9)/253431)*params[0]*delta + params[4]; //fieldoff
}

double vetopulsefitfunc::FitFunc_1st(double *x, double *params)
{
  double result = Eval(x, params);
  double delta = fMCFunc->GetBinWidth(1);
  double rate   = params[npar];
  double window = params[8];  
  double stRate = pow(rate,2)*window*exp(-rate*window);
  return result*(stRate)*params[0]*delta;
}

double vetopulsefitfunc::FitFunc_2nd(double *x, double *params)
{
  double result = Eval(x, params);
  double delta = fMCFunc->GetBinWidth(1);
  double rate   = params[npar];
  double window = params[8];
  double ndRate = pow(rate,3)*pow(window,2)*0.5*exp(-rate*window);
  return result*(ndRate)*params[0]*delta;
}
/*
double vetopulsefitfunc::SumFuncs(double *x,double *params)
{
  double sum = 0;
  for(size_t i=0; i<fitplot.size(); i++)
    {
      sum += fitplot.at(i)->EvalPar(x,params);
    }
  return sum;
}
*/






/*
double vetopulsefitfunc::C14Fit (double* x, double* params){
  int nBins = EnergyMC[0]->GetNbinsX();
  double delta = EnergyMC[0]->GetBinWidth(1);
  int choice = 1;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[0]->GetBinCenter(i);
      spectrum = EnergyMC[0]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double C14Rate = params[9];
  //  return result*(C14Rate)*delta;
  //  return result*(C14Rate*(1./6103.64)*1.e-7*177353.)*params[0]*delta;

  return result*(C14Rate*(2100.e-9)/449177)*params[0]*delta; //fieldon
  //  return result*(C14Rate*(2100.e-9)/253431)*params[0]*delta; //fieldoff
}

double vetopulsefitfunc::C14Fit_1st (double* x, double* params){
  int nBins = C14_EnergyMC_1st->GetNbinsX();
  double delta = C14_EnergyMC_1st->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  double window = params[8];
  for (int i=1; i<=nBins; i++)
    {
      energy = C14_EnergyMC_1st->GetBinCenter(i);
      spectrum = C14_EnergyMC_1st->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double stRate = pow(params[9],2)*window*exp(-params[9]*window);
  //  return result*(stRate)*delta;
  return result*(stRate)*params[0]*delta;
}

double vetopulsefitfunc::C14Fit_2nd (double* x, double* params){
  int nBins = C14_EnergyMC_2nd->GetNbinsX();
  double delta = C14_EnergyMC_2nd->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  double window = params[8];
  for (int i=1; i<=nBins; i++)
    {
      energy = C14_EnergyMC_2nd->GetBinCenter(i);
      spectrum = C14_EnergyMC_2nd->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double ndRate = pow(params[9],2)*window*exp(-params[9]*window);
  //  return result*(ndRate)*delta;
  return result*(ndRate)*params[0]*delta;
}

double vetopulsefitfunc::Co60Fit (double* x, double* params){
  int nBins = EnergyMC[1]->GetNbinsX();
  double delta = EnergyMC[1]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[1]->GetBinCenter(i);
      spectrum = EnergyMC[1]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double Co60Rate = params[10];
  //  return result*(Co60Rate)*delta;
  return result*(Co60Rate)*params[0]*delta;
}

double vetopulsefitfunc::Co57Fit (double* x, double* params){
  int nBins = EnergyMC[2]->GetNbinsX();
  double delta = EnergyMC[2]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[2]->GetBinCenter(i);
      spectrum = EnergyMC[2]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double Co57Rate = params[11];
  //  return result*(Co57Rate)*delta;
  return result*(Co57Rate)*params[0]*delta;
}

double vetopulsefitfunc::K40Fit (double* x, double* params){
  int nBins = EnergyMC[3]->GetNbinsX();
  double delta = EnergyMC[3]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[3]->GetBinCenter(i);
      spectrum = EnergyMC[3]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double K40Rate = params[12];
  //  return result*(K40Rate)*delta;
  return result*(K40Rate)*params[0]*delta;
}

double vetopulsefitfunc::Th232Fit (double* x, double* params){
  int nBins = EnergyMC[4]->GetNbinsX();
  double delta = EnergyMC[4]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[4]->GetBinCenter(i);
      spectrum = EnergyMC[4]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double Th232Rate = params[13];
  //  return result*(Th232Rate)*delta;
  return result*(Th232Rate)*params[0]*delta;
}
double vetopulsefitfunc::Th232LowerFit (double* x, double* params){
  int nBins = EnergyMC[5]->GetNbinsX();
  double delta = EnergyMC[5]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[5]->GetBinCenter(i);
      spectrum = EnergyMC[5]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double Th232LowerRate = params[14];
  return result*(Th232LowerRate)*params[0]*delta;
}

double vetopulsefitfunc::U235Fit (double* x, double* params){
  int nBins = EnergyMC[6]->GetNbinsX();
  double delta = EnergyMC[6]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[6]->GetBinCenter(i);
      spectrum = EnergyMC[6]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }

  double U235Rate = params[17]/21.5;
  return result*(U235Rate)*params[0]*delta;
}

double vetopulsefitfunc::U235LowerFit (double* x, double* params){
  int nBins = EnergyMC[7]->GetNbinsX();
  double delta = EnergyMC[7]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[7]->GetBinCenter(i);
      spectrum = EnergyMC[7]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double U235LowerRate = params[16];
  return result*(U235LowerRate)*params[0]*delta;
}

double vetopulsefitfunc::U238Fit (double* x, double* params){
  int nBins = EnergyMC[8]->GetNbinsX();
  double delta = EnergyMC[8]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[8]->GetBinCenter(i);
      spectrum = EnergyMC[8]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double U238UpperRate = params[17];
  //  return result*(U238UpperRate)*delta;
  return result*(U238UpperRate)*params[0]*delta;
}

double vetopulsefitfunc::U238LowerFit (double* x, double* params){
  int nBins = EnergyMC[9]->GetNbinsX();
  double delta = EnergyMC[9]->GetBinWidth(1);
  int choice = 0;
  double result = 0;
  double energy, spectrum;
  double q = x[0];
  for (int i=1; i<=nBins; i++)
    {
      energy = EnergyMC[9]->GetBinCenter(i);
      spectrum = EnergyMC[9]->GetBinContent(i);
      result += spectrum*response_function(q,energy,params,choice);
    }
  double U238LowerRate = params[18];
  //  return result*(U238LowerRate)*delta;
  return result*(U238LowerRate)*params[0]*delta;
}

double vetopulsefitfunc::TotalFit(Double_t *x, Double_t *params)
{
  double C14_value = C14Fit(x,params);
  //double C14_1st_value = C14Fit_1st(x,params);
  //  double C14_2nd_value = C14Fit_2nd(x,params);
  double Co60_value = Co60Fit(x,params);
  double Co57_value = Co57Fit(x,params);
  double K40_value = K40Fit(x,params);
  double Th232_value = Th232Fit(x,params);
  double Th232Lower_value = Th232LowerFit(x,params);
  double U235_value = U235Fit(x,params);
  double U235Lower_value = U235LowerFit(x,params);
  double U238_value = U238Fit(x,params);
  double U238Lower_value = U238LowerFit(x,params);

  return C14_value+Co60_value+Co57_value+K40_value+Th232_value+U238_value+U238Lower_value+U235_value+Th232Lower_value+U235Lower_value+params[4];

}

struct FitFunc
{
  FitFunc(TH1F* f, int k): fMCFunc(f), npar(k) {}

  double Eval(double *x, double *params){
      int nBins = fMCFunc->GetNbinsX();
      double delta = fMCFunc->GetBinWidth(1);
      int choice = 0;
      double result = 0;
      double energy, spectrum;
      double q = x[0];
      for (int i=1; i<=nBins; i++)
	{
	  energy   = fMCFunc->GetBinCenter(i);
	  spectrum = fMCFunc->GetBinContent(i);
	  result  += spectrum*response_function(q,energy,params,choice);
	}
      double rate = params[npar];
      return result*(rate)*params[0]*delta;    
  }
  
  TH1F* fMCFunc;
  int   npar;

};

*/
