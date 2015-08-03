#ifndef _vetopulsereal_H
#define _vetopulsereal_H

#include "TH1F.h"
#include "TH2F.h"
#include "TChain.h"
#include "TRint.h"
#include "TNtuple.h"
#include "TProof.h"
#include "TSelector.h"

#include "vetopulsebasic.hh"
#include "./odselector/odselector.h"
#include "./DSTtreeSelector/DSTtreeSelector.h"
#include "./SLADDSTSelector/SLADDSTSelector.h"

using namespace std;

class vetopulsebasic;
class DSTtreeSelector;
class odselector;
class SLADDSTSelector;

class vetopulsereal:public vetopulsebasic
{
public :
  vetopulsereal(){}

  virtual ~vetopulsereal()
  {}

  void DST_Readdatafile(TChain*, int, int, bool);
  void DST_Process(TChain*,TString);
  void DST_Selector(int,int,int);

  void OD_Readdatafile(TChain*, int, int);
  void OD_Process(TChain*,TString);
  void OD_Selector(int,int);

  void SLADDST_Readdatafile(TChain*, int, int);
  void SLADDST_Process(TChain*,TString);
  void SLADDST_Selector(int,int);
  
  bool multicut(float,float,float);

  ClassDef(vetopulsereal,0);

};

#endif /* _vetopulsereal_H */
