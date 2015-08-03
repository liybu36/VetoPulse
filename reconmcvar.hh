#ifndef reconmcvar_hh
#define reconmcvar_hh

#include <TObject.h>
#include "TROOT.h"
#include "TTree.h"
#include "TObject.h"
#include <string>
#include <vector>
#include <cstdlib>

using namespace std;

struct reconmcvar
{
  reconmcvar():
    pdg(0), event_broken(0), et(0), ex(0), ey(0), ez(0), edep(0),
    eqch(0),quenchingfactor(0),volume(0)
  {}

  virtual ~reconmcvar(){}
  
  Int_t           pdg;
  Bool_t          event_broken;
  vector<TString> *volume;
  vector<double>  *et;
  vector<double>  *ex;
  vector<double>  *ey;
  vector<double>  *ez;
  vector<double>  *edep;
  vector<double>  *eqch;
  vector<double>  *quenchingfactor;

};

#endif
