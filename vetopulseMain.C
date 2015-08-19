/*
  How to run it:
  set the isotopes in vetopulsebasic.cc file 
  make
  ./vetopulseMain cfg/testcfg.txt
 */
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

#include "vetopulsemc.hh"
#include "vetopulsereal.hh"
#include "vetopulsefitfunc.hh"
#include "vetopulsebasic.hh"
#include "vetopulsetotalfit.hh"
#include "vetopulseslad.hh"
#include "reconmcvar.hh"
#include "vetoreadcfg.hh"
#include "./odselector/odselector.h"
#include "./DSTtreeSelector/DSTtreeSelector.h"
#include "./SLADDSTSelector/SLADDSTSelector.h"

using namespace std;

TRint* theApp;

void vetopulseMain(char *cfile)
{
  Config config;
  config.loadConfig(cfile);
  cout<<config.DST<<" "<<config.OD<<" "<<config.sladdst<<endl;
    
  //  vetopulsebasic *basic = new vetopulsebasic;  
  //  basic->SetIsotopes();

  if(config.mc_tag){
    vetopulsemc *mc = new vetopulsemc();
    mc->SetMCinputdir(config.mcindir);
    mc->SetMCFile(config.mcfile);
    mc->SetBasics();
    
    mc->MCData_Analysis(config.mcu238);
    delete mc;
  }
  if(config.real_tag){
    vetopulsereal *real = new vetopulsereal();
    real->SetRealinputdir(config.realindir);
    real->SetRealFile(config.realfile);
    //    real->SetIsotopes();
   
    if(config.DST)
      real->DST_Selector(config.startrun,config.endrun,config.fieldon);
    if(config.OD)
      real->OD_Selector(config.startrun,config.endrun);
    // if(config.sladdst)
    // real->SLADDST_Selector(config.startrun,config.endrun); 
    delete real;
  }
  if(config.real_tag && config.sladdst)
    {
      vetopulseslad *slad = new vetopulseslad();
      slad->SetRealinputdir(config.realindir);
      slad->SetRealoutputdir(config.realoutputdir);
      slad->SetRealFile(config.realfile);
      slad->SLADProcess(config.startrun,config.endrun);

      delete slad;
    }
  if(config.totalfit_tag){
    vetopulsetotalfit *totalfit = new vetopulsetotalfit();
    totalfit->SetMCinputdir(config.mcindir);
    totalfit->SetMCFile(config.mcfile);
    totalfit->SetRealinputdir(config.realindir);
    totalfit->SetRealFile(config.realfile);
    totalfit->Setoutputdir(config.outdir);
    totalfit->SetOutFile(config.outfile);  
    totalfit->SetBasics();

    totalfit->TotalFit(config.startfit);
    delete totalfit;
  }
  

}

#ifndef __CINT__
int main(int argc, char** argv){
  theApp = new TRint("App",&argc,argv,NULL,0);
  char *cfile;
  if(theApp->Argc() == 2)
    cfile = theApp->Argv(1);
  else{
    cout<<"Usage: ./reconMain cfg.txt"<<endl;
    return 0;    
  }
  vetopulseMain(cfile);
  cout<<"!!! Application Processed Successfully "<<endl;
}
#endif
