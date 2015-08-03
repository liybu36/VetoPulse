#ifndef _vetopulsemc_H
#define _vetopulsemc_H
#include <map>
#include <vector>
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

#include "reconmcvar.hh"
#include "vetopulsebasic.hh"

using namespace std;

class vetopulsebasic;
struct reconmcvar;

class vetopulsemc:public vetopulsebasic
{
public :
  vetopulsemc():
    Histlist(0),MCValue_ntuple(0)
  {}
  
  virtual ~vetopulsemc()
  {}
  
  //  reconmcvar recon;
  vector<reconmcvar*> e; 

  void Recon_Readdatafile(TChain*,int,int,string,string);
  void MCData_Analysis(bool);  
  void Analyze_MCData(int);
  void Analyze_MCCoinData(int);
  void Analyze_U238MCCoinData(int,int,bool);
  void PDG_Analysis(int);
  void GetFraction(int,int);
  void BookHistograms();
  void Load_ReconTree(TChain*, reconmcvar*);
  void Fill_U_Pos();
  bool Check_U_Pos(int);

  void Calculate_Convolution(TH1F*, TH1F*, TH1F*);
  void Calculate_Norm(TH1F*, vector<double> &);

private:
  vector<TH1F*> EnergyMC;
  vector<TH1F*> EdepMC;
  vector<TH1F*> TPC_EnergyMC;
  vector<TH1F*> TPC_EdepMC;
  vector<TH2F*> TPC_Veto_EnergyMC;
  TNtuple* MCValue_ntuple;
  TObjArray Histlist;
  vector<int> U_Pos;

  ClassDef(vetopulsemc,0);
};

#endif /* _vetopulsemc_H */
