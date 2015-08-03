#include "vetopulsemc.hh"
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

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
#include "TCut.h"
#include "TLegend.h"
#include "TPad.h"
#include "TNtuple.h"
#include "TString.h"

using namespace std;

void vetopulsemc::Calculate_Convolution(TH1F *energy_spectrum1, TH1F *energy_spectrum2, TH1F *convolution)
{
  int np1 = energy_spectrum1->GetNbinsX();
  int np2 = energy_spectrum2->GetNbinsX();
  int nconv = np1 + np2;
  double binsize;
  for (int n = 1; n <= nconv; n++) {
    double sum = 0;
    for (int m = 1; m <= n; m++) {
      if (m > np1  || (n - m) > np2 ) continue;
      sum += (energy_spectrum1->GetBinContent(m)) * (energy_spectrum2->GetBinContent(n-m+1));
    }
    convolution->SetBinContent(n,sum);
  }
  double integral = convolution->Integral();
  cout<<convolution->GetName()<<"\t integral= "<<integral<<endl;
  convolution->Scale(1./integral);
}

void vetopulsemc::Calculate_Norm (TH1F *energy_spectrum,  vector<double> &convolution_spectrum)
{
  int np = energy_spectrum->GetNbinsX();
  double integral = energy_spectrum->Integral();
  cout<<"integral= "<<integral<<endl;
  energy_spectrum->Scale(1./integral);
  for(int i=1; i<=np; i++)
    convolution_spectrum.push_back(energy_spectrum->GetBinContent(i));
}

bool Compare_Time(double a, double b)
{
  return a<b;
}

void vetopulsemc::Load_ReconTree(TChain *recontree, reconmcvar* f)
{
  recontree->SetBranchAddress("volume",&f->volume);
  recontree->SetBranchAddress("et",&f->et);
  recontree->SetBranchAddress("ex",&f->ex);
  recontree->SetBranchAddress("ey",&f->ey);
  recontree->SetBranchAddress("ez",&f->ez);
  recontree->SetBranchAddress("edep",&f->edep);
  recontree->SetBranchAddress("eqch",&f->eqch);
  recontree->SetBranchAddress("quenchingfactor",&f->quenchingfactor);  
  recontree->SetBranchAddress("event_pdg", &f->pdg);
  recontree->SetBranchAddress("event_broken", &f->event_broken);

}

void vetopulsemc::GetFraction(int k, int NEntries)
{
  Int_t NBins = EnergyMC.at(k)->GetNbinsX();
  Double_t integralsum = EnergyMC.at(k)->Integral();
  cout<<"NBins= "<<NBins<<"\t integralsum="<<integralsum<<endl;
  EnergyMC.at(k)->Scale(1./Normalization.at(k));
  fraction.at(k) = 1.0*integralsum/NEntries;  
}

void vetopulsemc::Recon_Readdatafile(TChain *t,int start,int end,string middle,string last)
{
  string dirname=GetMCinputdir();
  for(int i=start; i<=end; i++)
    {
      TString filename;
      if(i==0) filename.Form("%s%s",middle.c_str(),last.c_str());
      else
	filename.Form("%s_v%d%s",middle.c_str(),i,last.c_str());
      filename.Prepend(dirname.c_str());
      if(Verifydatafile(filename))
	t->Add(filename);
    }
}

void vetopulsemc::Analyze_MCData(int k)
{
  TChain *recontree = new TChain("Recon");
  string middle = "out"+Source.at(k);
  string last = "_clustered.root";
  Recon_Readdatafile(recontree,4,4,middle,last);

  Int_t NEntries = recontree->GetEntries();
  cout<<"NEntries= "<<NEntries<<endl;  
  Load_ReconTree(recontree,e.at(k));

  for(Int_t i=0; i<NEntries; i++)
    { recontree->GetEntry(i);
      for(Int_t j=0; j<e.at(k)->et->size(); j++)
	{
	  if(e.at(k)->volume->at(j)=="p_scint")
	    {
	      EnergyMC.at(k)->Fill(e.at(k)->eqch->at(j));
	      EdepMC.at(k)->Fill(e.at(k)->edep->at(j));
	    }
	}
    }
  GetFraction(k,NEntries);
}

void vetopulsemc::Analyze_MCCoinData(int k)
{
  string middle = "out"+Source.at(k);
  string last="_clustered.root";
  TChain *recontree = new TChain("Recon");
  Recon_Readdatafile(recontree,4,4,middle,last);

  Int_t NEntries = recontree->GetEntries();
  cout<<"NEntries= "<<NEntries<<endl;
  Load_ReconTree(recontree,e.at(k));

  double tpc_low_threshold=0.5;
  double prompt_time = -50.; //ns
  double delay_time = 50.; //ns

  for(Int_t i=0; i<NEntries; i++)
    { recontree->GetEntry(i);
      double tpc_total=0;
      vector<double> tpc_trigger_time;
      for(Int_t j=0; j<e.at(k)->et->size(); j++)
	{
	  if(e.at(k)->volume->at(j)=="p_active" && e.at(k)->eqch->at(j)>tpc_low_threshold)
	    {
	      tpc_trigger_time.push_back(e.at(k)->et->at(j)*1.e+9); //ns
	      TPC_EnergyMC.at(k)->Fill(e.at(k)->eqch->at(j));
	      TPC_EdepMC.at(k)->Fill(e.at(k)->edep->at(j));
	      tpc_total += e.at(k)->eqch->at(j);
	    }
	}
      if(tpc_trigger_time.size()>1)
	std::sort(tpc_trigger_time.begin(),tpc_trigger_time.end(),Compare_Time);
      if(tpc_trigger_time.size())
	{
	  for(Int_t j=0; j<e.at(k)->et->size(); j++)
	    {
	      if(e.at(k)->volume->at(j)=="p_scint")
		{ double gps = e.at(k)->et->at(j)*1.e+9 - tpc_trigger_time.front();
		  // if(gps>prompt_time && gps<delay_time)
		  if(gps>-2000. && gps<100.)
		    {
		      EnergyMC.at(k)->Fill(e.at(k)->eqch->at(j));
		      EdepMC.at(k)->Fill(e.at(k)->edep->at(j));
		      TPC_Veto_EnergyMC.at(k)->Fill(tpc_total,e.at(k)->eqch->at(j));
		    }
		}
	    }
	}
    }
  EnergyMC.at(k)->Rebin(2);
  GetFraction(k,NEntries);
}

void vetopulsemc::Analyze_U238MCCoinData(int k, int m, bool coin)
{
  string middle = "out"+Source.at(k);
  string last="_clustered.root";
  TChain *recontree = new TChain("Recon");
  Recon_Readdatafile(recontree,4,4,middle,last);

  Int_t NEntries = recontree->GetEntries();
  cout<<"NEntries= "<<NEntries<<endl;
  Load_ReconTree(recontree,e.at(k));  

  double tpc_low_threshold=0.5;
  double prompt_time = -50.; //ns
  double delay_time = 50.; //ns

  for(Int_t i=0; i<NEntries; i++)
    { recontree->GetEntry(i);
      vector<double> tpc_trigger_time;
      for(Int_t j=0; j<e.at(k)->et->size(); j++)
	{
	  if(e.at(k)->volume->at(j)=="p_active" && e.at(k)->eqch->at(j)>tpc_low_threshold)
	    {
	      tpc_trigger_time.push_back(e.at(k)->et->at(j)*1.e+9);
	      if(!e.at(k)->event_broken)
		{
		  TPC_EnergyMC.at(k)->Fill(e.at(k)->eqch->at(j));
		  TPC_EdepMC.at(k)->Fill(e.at(k)->edep->at(j));
		}
	      else
		{
		  TPC_EnergyMC.at(m)->Fill(e.at(k)->eqch->at(j));
		  TPC_EdepMC.at(m)->Fill(e.at(k)->edep->at(j));
		}
	    }
	}
      if(coin){
      if(tpc_trigger_time.size()>1)
	std::sort(tpc_trigger_time.begin(),tpc_trigger_time.end(),Compare_Time);
      if(tpc_trigger_time.size())
	{
	  for(Int_t j=0; j<e.at(k)->et->size(); j++)
	    {
	      if(e.at(k)->volume->at(j)=="p_scint")
		{ double gps = e.at(k)->et->at(j)*1.e+9 - tpc_trigger_time.front();
		  if(gps>-2000. && gps<100.)
		    //  if(gps>prompt_time && gps<delay_time)
		    {
		      if(!e.at(k)->event_broken)
			{
			  EnergyMC.at(k)->Fill(e.at(k)->eqch->at(j));
			  EdepMC.at(k)->Fill(e.at(k)->edep->at(j));
			}
		      else{
			EnergyMC.at(m)->Fill(e.at(k)->eqch->at(j));
			EdepMC.at(m)->Fill(e.at(k)->edep->at(j));
		      }  }
		} } }
      }else{
      for(Int_t j=0; j<e.at(k)->et->size(); j++)
	{
	  if(e.at(k)->volume->at(j)=="p_scint")
	    {
	      if(!e.at(k)->event_broken){
		EnergyMC.at(k)->Fill(e.at(k)->eqch->at(j));
		EdepMC.at(k)->Fill(e.at(k)->edep->at(j));
	      }
	      else{
		EnergyMC.at(m)->Fill(e.at(k)->eqch->at(j));
		EdepMC.at(m)->Fill(e.at(k)->edep->at(j));
	      }
	    }
	}
      }
    }
  EnergyMC.at(k)->Rebin(2);
  GetFraction(k,NEntries);

  EnergyMC.at(m)->Rebin(2);
  GetFraction(m,NEntries);  
}

void vetopulsemc::PDG_Analysis(int k)
{
  string middle = "out"+Source.at(k);
  string last=".root";
  TChain *dstree = new TChain("dstree");
  Recon_Readdatafile(dstree,0,6,middle,last);

  Int_t NEntries = dstree->GetEntries();
  cout<<"NEntries= "<<NEntries<<endl;

  Int_t pdg;
  dstree->SetBranchAddress("pdg", &pdg);

  std::map<Int_t,Int_t> pdg_counts;
  for(int i=0; i<NEntries; i++)
    {
      dstree->GetEntry(i);
      if(pdg_counts.count(pdg)==0)
	pdg_counts.insert( std::pair<Int_t,Int_t>(pdg,1) );
      else{
	++(pdg_counts.find(pdg)->second);
      }
    }
  for (std::map<Int_t, Int_t>::iterator it=pdg_counts.begin(); it!=pdg_counts.end(); ++it)
    cout << it->first << " => " << it->second << '\n';
  Normalization.at(k) = pdg_counts.find(PDG.at(k))->second;
  cout<<"Normalization= "<<Normalization.at(k)<<endl;

  if(Check_U_Pos(k))
    {
      Normalization.at(k+1) = pdg_counts.find(PDG.at(k+1))->second;
      cout<<"Normalization= "<<Normalization.at(k+1)<<endl;
    }
}

void vetopulsemc::BookHistograms()
{
  string temp = "_EnergyMC";
  for(size_t i=0; i<Source.size(); i++)
    {
      e.push_back(new reconmcvar());
    }  

  for(size_t i=0; i<Source_Pos.size(); i++)
    {
      int k = Source_Pos.at(i);
      //     e.push_back(new reconmcvar());
      EnergyMC.push_back(new TH1F(Form("%s%s",Source.at(k).c_str(),temp.c_str()),
				  Form("%s MC Energy",Source.at(k).c_str()),NBins.at(k),0,NBins.at(k)*1.0));
      EnergyMC.back()->GetXaxis()->SetTitle("Quenching Energy [keVee]");
      Histlist.Add(EnergyMC.back());

      EdepMC.push_back(new TH1F(Form("%s_EdepMC",Source.at(k).c_str()),"; Energy [keV]",NBins.at(k),0,NBins.at(k)*1.0));
      Histlist.Add(EdepMC.back());

      TPC_EnergyMC.push_back(new TH1F(Form("%s_TPC%s",Source.at(k).c_str(),temp.c_str()),
				      Form("%s TPC MC Energy",Source.at(k).c_str()),TBins.at(k),0,TBins.at(k)*1.0));
      Histlist.Add(TPC_EnergyMC.back());
      
      TPC_EdepMC.push_back(new TH1F(Form("%s_TPC_EdepMC",Source.at(k).c_str()),"; Energy [keV]",TBins.at(k),0,TBins.at(k)*1.0));
      Histlist.Add(TPC_EdepMC.back());

      TPC_Veto_EnergyMC.push_back(new TH2F(Form("%s_TPC_Veto%s",Source.at(k).c_str(),temp.c_str()),
					   Form("%s TPC Veto MC Energy",Source.at(k).c_str()),NBins.at(k),0,NBins.at(k)*1.0,NBins.at(k),0,NBins.at(k)*1.0));
      TPC_Veto_EnergyMC.back()->GetXaxis()->SetTitle("TPC Energy [keV]");
      TPC_Veto_EnergyMC.back()->GetYaxis()->SetTitle("Veto Energy [keV]");
      Histlist.Add(TPC_Veto_EnergyMC.back());
    }
  MCValue_ntuple = new TNtuple("MCValue_ntuple","MC Stats Value","fraction:Normalization");
  Histlist.Add(MCValue_ntuple);
}

void vetopulsemc::Fill_U_Pos()
{
  U_Pos.push_back(Search_Isotope("U235"));
  U_Pos.push_back(Search_Isotope("U238"));
}

bool vetopulsemc::Check_U_Pos(int k)
{
  for(size_t i=0; i<U_Pos.size(); ++i)
    if(k==U_Pos.at(i))
      return true;
}

void vetopulsemc::MCData_Analysis(bool coin)
{
  Fill_U_Pos();
  BookHistograms();
  for(size_t i=0; i<Source_Pos.size(); i++)
    {
      int k = Source_Pos.at(i);
      if(k==0)
	//if(k==0 || k==1 || k==3)
	{
	  PDG_Analysis(k);
	  Analyze_MCData(k);
	}
      else if(Check_U_Pos(k)) //U235 and U238 Chain
	{
	  PDG_Analysis(k);
	  Analyze_U238MCCoinData(k,k+1,coin);
	  ++i;
	}
      else{
	PDG_Analysis(k);
	Analyze_MCCoinData(k);     
      }
    }
  for(size_t i=0; i<Source_Pos.size(); i++)
    MCValue_ntuple->Fill(fraction.at(i),Normalization.at(i));

  string output = GetMCinputdir() + GetMCFile();
  cout<<output<<endl;
  TFile outfile(output.c_str(), "RECREATE");
  Histlist.Write();
  outfile.Write();
  outfile.Close();
  cout<<"Successfully save MC Data."<<endl;
}
