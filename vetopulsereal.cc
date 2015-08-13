#include "vetopulsereal.hh"
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
#include "TSelector.h"

//#include "./odselector/odselector.h"
//#include "./DSTtreeSelector/DSTtreeSelector.h"

using namespace std;

bool vetopulsereal::multicut(float height,float multiplicity, float charge){
  return height/multiplicity < (2.563e7 + TMath::Sqrt(1.574e14+1.390e12*(charge-14.40)*(charge-14.40)));
}

//Run the DST File
void vetopulsereal::DST_Readdatafile(TChain *t, int start, int end, bool fieldon)
{
  if(fieldon)
    {
      //run 11856 to 11928
      //      string dirname="/darkside/data/UAr_DSTs/";
      string dirname="/darkside/users/hqian/pulse_splitter/UAr_masadata/";
      //string dirname=GetRealinputdir();
      string middle="DST_Run";
      string last=".root";
      for(int i=start; i<=end; i++)
	{
	  TString filename = Form("%s%06d%s",middle.c_str(),i,last.c_str());
	  filename.Prepend(dirname.c_str());
	  if(Verifydatafile(filename))
	    t->Add(filename);
	}
    }
  else{
    string dirname="/darkside/data/UAr_DSTs/nullfield/";
    TString filename = dirname + "DST_Run011764_11822.root";
    if(Verifydatafile(filename))
      t->Add(filename);    
  }
}

void vetopulsereal::DST_Process(TChain* chain, TString label)
{
    TString option = label;
#define DST_TProof
#ifdef  DST_TProof
    chain->SetProof();
    TProof* pr = TProof::Open("workers=2");
    pr->SetParameter("PROOF_Packetizer","TPacketizer");
    pr->SetParameter("PROOF_MaxSlavesPerNode",8);
    chain->Process("./DSTtreeSelector/DSTtreeSelector.C+",option.Data(),-1,0);
#else
    DSTtreeSelector *selector = new DSTtreeSelector();
    selector->SetTRint(theApp);
    chain->Process(selector,option.Data());
    gSystem->Exit(0);    
#endif
}

void vetopulsereal::DST_Selector(int start, int end,int fieldon)
{
  TChain *DSTtree =new TChain("DSTtree");
  DST_Readdatafile(DSTtree,start,end,fieldon);
  Int_t nEntries = DSTtree->GetEntries();
  string label = GetRealFile();
  cout<<"nEntries= "<<nEntries<<" file:"<<GetRealFile()<<endl;
  DST_Process(DSTtree,label);
}


//Use Selector for ODRun000000.root files
void vetopulsereal::OD_Readdatafile(TChain *t, int start, int end)
{
  string dirname=GetRealinputdir();
  string middle="ODRun";
  string last=".root";
  TString filename;

  for(int i=start; i<=end; i++)
    {
      filename.Form("%s%s%06d%s",dirname.c_str(),middle.c_str(),i,last.c_str());
      if(Verifydatafile(filename))
	t->Add(filename);
    }
}

void vetopulsereal::OD_Process(TChain* chain, TString label){
  TString option = label;
#define OD_TProof
#ifdef  OD_TProof
  chain->SetProof();
  TProof* pr = TProof::Open("workers=2");
  pr->SetParameter("PROOF_Packetizer","TPacketizer");
  pr->SetParameter("PROOF_MaxSlavesPerNode",8);
  chain->Process("./odselector/odselector.C+",option.Data(),-1,0);
#else
  odselector *selector = new odselector();
  selector->SetTRint(theApp);
  chain->Process(selector,option.Data());
  gSystem->Exit(0);
#endif
}

void vetopulsereal::OD_Selector(int start, int end)
{
  TChain *odtree =new TChain("odtree");
  OD_Readdatafile(odtree,start,end);
  Int_t nEntries = odtree->GetEntries();
  cout<<"nEntries= "<<nEntries<<endl;
  string label = GetRealinputdir() + GetRealFile();
  OD_Process(odtree,label);
}


//Use Selector for SLADDST_Run000000.root files
void vetopulsereal::SLADDST_Readdatafile(TChain *t, int start, int end)
{
  string dirname=GetRealinputdir();
  string middle="SLADDST_Run";
  string last=".root";
  TString filename;

  for(int i=start; i<=end; i++)
    {
      filename.Form("%s%s%06d%s",dirname.c_str(),middle.c_str(),i,last.c_str());
      if(Verifydatafile(filename))
	t->Add(filename);
    }
}

void vetopulsereal::SLADDST_Process(TChain* chain, TString label){
  TString option = label;
#define SLADDST_TProof
#ifdef  SLADDST_TProof
  chain->SetProof();
  TProof* pr = TProof::Open("workers=2");
  pr->SetParameter("PROOF_Packetizer","TPacketizer");
  pr->SetParameter("PROOF_MaxSlavesPerNode",8);
  chain->Process("./SLADDSTSelector/SLADDSTSelector.C+",option.Data(),-1,0);
#else
  SLADDSTSelector *selector = new SLADDSTSelector();
  selector->SetTRint(theApp);
  chain->Process(selector,option.Data());
  gSystem->Exit(0);
#endif
}

void vetopulsereal::SLADDST_Selector(int start, int end)
{
  TChain *DSTtree =new TChain("DSTtree");
  SLADDST_Readdatafile(DSTtree,start,end);
  Int_t nEntries = DSTtree->GetEntries();
  cout<<"nEntries= "<<nEntries<<endl;
  string label = GetRealFile();
  SLADDST_Process(DSTtree,label);
}

