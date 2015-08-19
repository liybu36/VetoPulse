#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TRint.h"
#include "TProof.h"
#include "TCut.h"
#include "TLegend.h"
#include "vetopulsetotalfit.hh"

using namespace std;

void vetopulsetotalfit::Init()
{
  startrange = 0.;
  endrange   = 2000.;
  primary_npar = 9;  

}

double vetopulsetotalfit::Isotope_Sum(int isotopes)
{
  double activity_sum=0;
  cout<<isotopes<<"  "<<isotope_activity.size()<<endl;
  for(int i=0; i<isotopes; i++)
    {
      activity_sum += isotope_activity.at(i);       
    }
  cout<<"Calculating Isotope_Sum"<<endl;
  return activity_sum;
}

void vetopulsetotalfit::Isotope_Fraction(int isotopes)
{
  double sum = Isotope_Sum(isotopes);
  for(int i=0; i<isotopes; i++)
    {
      isotope_fraction.push_back(isotope_activity.at(i)*1.0/sum);
    }
  cout<<"Calculating Isotope_Fraction"<<endl;
}

bool vetopulsetotalfit::Isotope_Activity(int isotopes)
{
  MCValue_ntuple = (TNtuple*) MCFile->Get("MCValue_ntuple");

  float fraction;
  MCValue_ntuple->SetBranchAddress("fraction",&fraction);
  int nEntries = MCValue_ntuple->GetEntries();
  for(int i=0; i<nEntries; i++)
    {
      MCValue_ntuple->GetEntry(i);
      if(fraction>0)
	isotope_activity.push_back(fraction);
      else continue;
    }
  cout<<"Calculating Isotope_Activity "<<isotope_activity.size()<<endl;
  Isotope_Fraction(isotopes);  
  return true;
}

bool vetopulsetotalfit::Each_Isotope_Sum(int isotopes)
{
  if(!Isotope_Activity(isotopes))
    return false;
  for(int i=0; i<isotopes; i++)
    {
      integral_sum.push_back(isotope_fraction.at(i)*integralsum);
      cout<<"Source= "<<Source.at(Source_Pos.at(i))<<"\t isotope_fraction[i] "<<isotope_fraction.at(i)<<
	"\t integral_sum="<<integral_sum.at(i)<<endl;
    }
  cout<<"Calculating Each_Isotope_Sum"<<endl;
  return true;
}

double vetopulsetotalfit::GetIntegral(double startrange, double endrange)
{
  int Bins = FullSpectrum->GetNbinsX();
  int startbin = FullSpectrum->FindBin(startrange);
  int endbin = FullSpectrum->FindBin(endrange);
  double sum = FullSpectrum->Integral(startbin,endbin);
  cout<<"total integralsum= "<<sum<<endl;
  return sum;
}

double vetopulsetotalfit::GetLiveTime()
{
  TNtuple *LiveTime = (TNtuple*)RealFile->Get("livetime");
  float sum = 0;
  float live_time;
  LiveTime->SetBranchAddress("live_time",&live_time);
  for(int i=0; i<LiveTime->GetEntries(); i++)
    {
      LiveTime->GetEntry(i);
      sum += live_time;
    }
  cout<<"Live Time is "<<sum<<endl;
  return sum;
}

void vetopulsetotalfit::ScaleSpectrum(TH1F* f)
{
  int rebins = 5;
  f->Rebin(rebins);
  f->Sumw2();
  double livetime = GetLiveTime();
  /*
#ifdef fieldon
  double livetime = 449177; //6103.64; //[s]
#else
  double livetime = 253431; //[s]
#endif
  */
  f->Scale(1./livetime);
  f->SetTitle(""); 
}

bool vetopulsetotalfit::Load_MCPlots(int isotopes)
{
  MCDataFile = GetMCinputdir() + GetMCFile();
  if(!Verifydatafile(MCDataFile)){
    cout<<"!!!Error: Cannot Load_MCPlots"<<endl;
    return false;
  }
  MCFile = new TFile(MCDataFile.c_str());
  for(int i=0; i<isotopes; i++){
    mcplot.push_back((TH1F*)MCFile->Get(Form("%s_EnergyMC",Source.at(Source_Pos.at(i)).c_str())));
    MClist.Add(mcplot.at(i));
  }  
  return true;
}

bool vetopulsetotalfit::Load_RealPlots()
{
  string RealDataFile = GetRealinputdir() + GetRealFile();
  cout<<"Fitting: "<<endl;
  if(!Verifydatafile(RealDataFile)){
    cout<<"!!!Error: Cannot Load_RealPlots"<<endl;
    return false;
  }
  RealFile = new TFile(RealDataFile.c_str());
  FullSpectrum = (TH1F*) RealFile->Get("FullSpectrum");
  return true;
}

void vetopulsetotalfit::SetPrimaryVar()
{
  Double_t binw             = FullSpectrum->GetBinWidth(1);
  Double_t ly_mean          = 0.6058;  //0.63;
  Double_t spe_var          = 0.14;
  Double_t ly_var           = 0.0043;  // 0.00642;
  Double_t constant         = 0.;
  Double_t basel_mean       = 0.386;
  Double_t basel_var        = 100; //1.;
  Double_t threshold        = 0; //1;
  Double_t window           = 5.78e-7;   // 2e-7;             

  parvalues.push_back(binw);
  parvalues.push_back(ly_mean);
  parvalues.push_back(spe_var);
  parvalues.push_back(ly_var);
  parvalues.push_back(constant);
  parvalues.push_back(basel_mean);
  parvalues.push_back(basel_var);
  parvalues.push_back(threshold);
  parvalues.push_back(window);

  parnames.push_back("Bin Width");
  parnames.push_back("LY Mean[PE/keV]");
  parnames.push_back("Rel SPE Var");
  parnames.push_back("Rel LY Var");
  parnames.push_back("Constant");
  parnames.push_back("Baseline Mean[PE]");
  parnames.push_back("Baseline Var");
  parnames.push_back("Threshold");
  parnames.push_back("Window With[s]");
 
  paruplimits.push_back(20);
  paruplimits.push_back(1);
  paruplimits.push_back(1);
  paruplimits.push_back(1);
  paruplimits.push_back(10);
  paruplimits.push_back(50);
  paruplimits.push_back(100);
  paruplimits.push_back(100);
  paruplimits.push_back(1.e-6);

}

void vetopulsetotalfit::SetFitVar(int isotopes)
{
  for(int i=0; i<isotopes; i++)
    {
      parvalues.push_back(integral_sum.at(i));
      parnames.push_back(Source.at(Source_Pos.at(i))+" Rate[Bq]");
      paruplimits.push_back(integralsum*100);
    }  
  cout<<"Set FixVar"<<endl;
}

double vetopulsetotalfit::SumFuncs(double *x, double *params)
{
  double sum = 0;
  for(size_t i=0; i<fitplot.size(); i++)
    {
      sum += fitplot.at(i)->EvalPar(x,params);
    }
  return sum+params[4];
}

bool vetopulsetotalfit::BookFitFuncs(int isotopes, int parnums)
{ 
  for(int i=0; i<isotopes; i++){
    func.push_back(new vetopulsefitfunc(mcplot.at(i),i+primary_npar));
    string fitplot_name = "Fit_" + Source.at(Source_Pos.at(i));
    /*
      if(Source.at(Source_Pos.at(i))=="Tl208")
      fitplot.push_back(new TF1(fitplot_name.c_str(),func.at(i),&vetopulsefitfunc::FitFunc,1000.,endrange,parnums));
      else
    */
    fitplot.push_back(new TF1(fitplot_name.c_str(),func.at(i),&vetopulsefitfunc::FitFunc,startrange,endrange,parnums));
    Fitlist.Add(fitplot.back());
  }
  //  TF1* Fit_C14 = new TF1("Fit_C14",split,&vetopulsesplit::C14Fit, startrange, endrange, parnums,"vetopulsesplit","C14Fit");
  Fit_Total = new TF1("Fit_Total",this,&vetopulsetotalfit::SumFuncs,startrange,endrange,parnums);
  
  return true;
}

void vetopulsetotalfit::FillFitVars(int parnums)
{
  for(int i=0; i<parnums; i++)
    {
      Fit_Total->SetParName(i,parnames[i].c_str());
      if(i==0||i==2||i==3||i==4||i==7||i==8||i==5||i==parnums-4){
	Fit_Total->FixParameter(i,parvalues[i]);
	cout<<"Fixed parnames: "<<i<<" "<<parnames.at(i)<<endl;
      }else{
	  Fit_Total->SetParameter(i, parvalues[i]);
	  Fit_Total->SetParLimits(i,0,paruplimits[i]);
	}
    }  
}

bool vetopulsetotalfit::TotalFit(int startfit)
{
  Init();
  int isotope_number = static_cast<int> (Source_Pos.size());
  cout<<"Number of isotopes: "<<isotope_number<<endl;

  if(!Load_MCPlots(isotope_number) || !Load_RealPlots())
    return false;

  ScaleSpectrum(FullSpectrum);
  TLegend *leg = new TLegend(0.5,0.7,0.8,1.0);
  canv.push_back(new TCanvas("c2", "Full Spectrum with Fit", 1000, 400));
  canv.back()->SetLogy();
  canv.back()->cd();
  FullSpectrum->Draw();

  integralsum = GetIntegral(startrange,endrange);  
  const Int_t parnums = isotope_number + primary_npar;

  if(!Each_Isotope_Sum(isotope_number))
    return false;
  SetPrimaryVar();
  SetFitVar(isotope_number);

  if(!BookFitFuncs(isotope_number,parnums))
    return false;

  FillFitVars(parnums);
  Fit_Total->SetNpx(2000);
  Fit_Total->SetLineColor(Colors(12));
  if(startfit)
    {
      //FullSpectrum->Fit(Fit_Total,"RVWLM");
      FullSpectrum->Fit(Fit_Total,"RVEM");
    }

  Double_t new_par[parnums];
  Fit_Total->GetParameters(new_par);
  for(Int_t i=0; i<Fitlist.GetEntries(); i++)
    {
      TF1* tempfit = dynamic_cast<TF1*>(Fitlist.At(i));
      tempfit->SetLineColor(Colors(i));
      tempfit->SetParameters(new_par);
      tempfit->Draw("SAME");
      leg->AddEntry(tempfit,tempfit->GetName(),"l");
      cout<<i<<"\t"<<tempfit->GetName()<<"  "<<Fitlist.GetSize()<<" "<<Fitlist.GetEntries()<<endl;
    }
  Fit_Total->Draw("same");
  leg->AddEntry(Fit_Total,"Total Fit","l");
  leg->Draw();
  
  Double_t chi2 = Fit_Total->GetChisquare();
  Int_t ndf = Fit_Total->GetNDF();
  cout<<"chi2/ndf= "<<chi2/(1.0*ndf)<<"\t Probability="<<Fit_Total->GetProb()<<endl;
  
  SaveHistograms();
  return true;
}

void vetopulsetotalfit::SaveHistograms()
{
  string output = Getoutputdir() + GetOutFile();
  TFile outfile(output.c_str(), "RECREATE");
  MClist.Write();
  Fitlist.Write();
  Fit_Total->Write();
  FullSpectrum->Write();
  for(size_t i=0; i<canv.size(); i++)
    canv.at(i)->Write();

  outfile.Write();
  outfile.Close();
}




  /*
  Double_t binw             = FullSpectrum->GetBinWidth(1);
  Double_t ly_mean          = 0.6058;  //0.63;
  Double_t spe_var          = 0.14;
  Double_t ly_var           = 0.0043;  // 0.00642;
  Double_t constant         = 0.;
  Double_t basel_mean       = 0.386;
  Double_t basel_var        = 100; //1.;
  Double_t threshold        = 0; //1;
  Double_t window           = 5.78e-7;   // 2e-7;             

  Double_t parvaluesamples[]={binw,ly_mean,spe_var,ly_var,constant,basel_mean,basel_var,threshold,window,
			      //integralsum,integralsum,integralsum,integralsum,integralsum,integralsum,integralsum};
			      integral_sum[0],integral_sum[1],integral_sum[2],integral_sum[3],
			      integral_sum[4],integral_sum[5],0,integral_sum[7],integral_sum[8]*10,integral_sum[9]*10};
  string parnamesamples[]={"Bin Width","LY Mean [PE/keV]","Rel SPE Var","Rel LY Var",
			 "Constant","Baseline Mean [p.e]","Baseline Var","Threshold","Window Width [s]",
			 "C14 Rate[Bq]","Co60 Rate[Bq]","Co57 Rate[Bq]","K40 Rate[Bq]",
			 "Th232Upper Rate[Bq]","Th232Lower Rate[Bq]","U235Upper Rate[Bq]","U235Lower Rate[Bq]","U238Upper Rate[Bq]","U238Lower Rate[Bq]"};
  // Double_t parlowlimits[parnums] = {0,0,0,0,-100,-10,-100,0,0,0,0,0,0,0,0};
  // Double_t paruplimits[]  = {10,1.,0.14,100,100,100,100,100,1e-6,integralsum*1,integralsum*2,integralsum*2,integralsum*2,integralsum*2,integralsum*2};
  
  Double_t parlowlimits[] = {0,0,0,0,0,0,0,0,0,
			     // 1,1,1,1,1,1,1};
			     0,0,0,0,0,0,0,0,0,0};
  Double_t paruplimits[]  = {20,1,1,1,10,50,100,100,1e-6,
			     integralsum*20,integralsum*1,integralsum*20,integralsum*20,
			     integralsum*100,integralsum*100,integralsum*100/21.5,integralsum*100,integralsum*100,integralsum*100};
    Fit_Total_name.erase(Fit_Total_name.rfind("+"));

  GetIsotopeActivity(int isotope_number);

  Double_t isotope_fraction[isotope_number];
  Isotope_Fraction(isotope_number,isotope_activity,isotope_fraction);
  Double_t integral_sum[isotope_number];  
  vector<int> linecolor = Colors();

  */

