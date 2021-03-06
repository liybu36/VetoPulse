//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun 25 14:31:18 2015 by ROOT version 5.34/25
// from TTree DSTtree/tree of selected events
// found on file: SLADDST_Run011138.root
//////////////////////////////////////////////////////////

#ifndef SLADDSTSelector_h
#define SLADDSTSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TNtuple.h>
#include "TProofOutputFile.h"
// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
#include <vector>
using namespace std;
// Fixed size dimensions of array or collections stored in the TTree if any.

class SLADDSTSelector : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   Int_t           runID;
   Int_t           tpc_eventID;
   Int_t           tpc_event_type;
   Double_t        tpc_gps_fine;
   Double_t        tpc_gps_coarse;
   Double_t        tpc_s1_start_time;
   Double_t        tpc_total_s1;
   Double_t        tpc_total_s1_f90;
   Double_t        tpc_total_s2;
   Double_t        tpc_total_s2_f90;
   Double_t        tpc_t_drift;
   Double_t        tpc_s1_late;
   Int_t           tpc_npulses;
   Double_t        tpc_timestamp;
   Double_t        tpc_live_time;
   Double_t        tpc_inhibit_time;
   Int_t           tpc_has_s3;
   Float_t         tpc_masa_x;
   Float_t         tpc_masa_y;
   Float_t         tpc_jason_x;
   Float_t         tpc_jaosn_y;
   Int_t           od_eventID;
   Double_t        od_gps_fine;
   Double_t        od_gps_coarse;
   Double_t        od_timestamp;
   Int_t           od_nclusters;
   Double_t        od_lsv_charge;
   Double_t        od_wt_charge;
   vector<double>  *od_cluster_charge;
   vector<double>  *od_cluster_start;
   vector<double>  *od_cluster_height;
   vector<double>  *od_cluster_multiplicity;
   vector<int>     *od_cluster_pass_multcut;
   vector<double>  *od_cluster_dtprompt;
   vector<int>     *od_roi_lsv_id;
   vector<float>   *od_roi_lsv_charge;
   vector<int>     *od_roi_lsv_max_multiplicity;
   vector<int>     *od_slider_lsv_id;
   vector<float>   *od_slider_lsv_charge;
   vector<float>   *od_slider_lsv_time;
   vector<int>     *od_slider_lsv_max_multiplicity;

   // List of branches
   TBranch        *b_runID;   //!
   TBranch        *b_tpc_eventID;   //!
   TBranch        *b_tpc_event_type;   //!
   TBranch        *b_tpc_gps_fine;   //!
   TBranch        *b_tpc_gps_coarse;   //!
   TBranch        *b_tpc_s1_start_time;   //!
   TBranch        *b_tpc_total_s1;   //!
   TBranch        *b_tpc_total_s1_f90;   //!
   TBranch        *b_tpc_total_s2;   //!
   TBranch        *b_tpc_total_s2_f90;   //!
   TBranch        *b_tpc_t_drift;   //!
   TBranch        *b_tpc_s1_late;   //!
   TBranch        *b_tpc_npulses;   //!
   TBranch        *b_tpc_timestamp;   //!
   TBranch        *b_tpc_live_time;   //!
   TBranch        *b_tpc_inhibit_time;   //!
   TBranch        *b_tpc_has_s3;   //!
   TBranch        *b_tpc_masa_x;   //!
   TBranch        *b_tpc_masa_y;   //!
   TBranch        *b_tpc_jason_x;   //!
   TBranch        *b_tpc_jason_y;   //!
   TBranch        *b_od_eventID;   //!
   TBranch        *b_od_gps_fine;   //!
   TBranch        *b_od_gps_coarse;   //!
   TBranch        *b_od_timestamp;   //!
   TBranch        *b_od_nclusters;   //!
   TBranch        *b_od_lsv_charge;   //!
   TBranch        *b_od_wt_charge;   //!
   TBranch        *b_od_cluster_charge;   //!
   TBranch        *b_od_cluster_start;   //!
   TBranch        *b_od_cluster_height;   //!
   TBranch        *b_od_cluster_multiplicity;   //!
   TBranch        *b_od_cluster_pass_multcut;   //!
   TBranch        *b_od_cluster_dtprompt;   //!
   TBranch        *b_od_roi_lsv_id;   //!
   TBranch        *b_od_roi_lsv_charge;   //!
   TBranch        *b_od_roi_lsv_max_multiplicity;   //!
   TBranch        *b_od_slider_lsv_id;   //!
   TBranch        *b_od_slider_lsv_charge;   //!
   TBranch        *b_od_slider_lsv_time;   //!
   TBranch        *b_od_slider_lsv_max_multiplicity;   //!

   double charge,width,height;
   
 SLADDSTSelector(TTree * /*tree*/ =0) : fChain(0),
     id(0),
     FullSpectrum(0),FullSpectrum_ntuple(0),C14Spectrum(0),Veto_TPC_hist(0),K40Spectrum(0),
     charge(0),height(0),width(0),Fulllate(0),Co60Spectrum(0),RandomSpectrum(0),
     F1Spectrum(0), F2Spectrum(0),lsvtree(0),fFile(0),fProofFile(0)
     { }
   virtual ~SLADDSTSelector() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

 private:
   void BookHistograms();
   void FillHistograms();
   TTree *lsvtree;
   TFile *fFile;
   TProofOutputFile *fProofFile;
   void SetOutputTree(TFile *);
   bool multicut(float,float,float);
   TH1F *id;
   TH1F *FullSpectrum, *Fulllate, *RandomSpectrum;
   TNtuple *FullSpectrum_ntuple;
   TH1F *C14Spectrum,*Co60Spectrum, *K40Spectrum;
   TH1F *F1Spectrum, *F2Spectrum;
   TH2F *Veto_TPC_hist;

   ClassDef(SLADDSTSelector,0);
};

#endif

#ifdef SLADDSTSelector_cxx
void SLADDSTSelector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   od_cluster_charge = 0;
   od_cluster_start = 0;
   od_cluster_height = 0;
   od_cluster_multiplicity = 0;
   od_cluster_pass_multcut = 0;
   od_cluster_dtprompt = 0;
   od_roi_lsv_id = 0;
   od_roi_lsv_charge = 0;
   od_roi_lsv_max_multiplicity = 0;
   od_slider_lsv_id = 0;
   od_slider_lsv_charge = 0;
   od_slider_lsv_time = 0;
   od_slider_lsv_max_multiplicity = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runID", &runID, &b_runID);
   fChain->SetBranchAddress("tpc_eventID", &tpc_eventID, &b_tpc_eventID);
   fChain->SetBranchAddress("tpc_event_type", &tpc_event_type, &b_tpc_event_type);
   fChain->SetBranchAddress("tpc_gps_fine", &tpc_gps_fine, &b_tpc_gps_fine);
   fChain->SetBranchAddress("tpc_gps_coarse", &tpc_gps_coarse, &b_tpc_gps_coarse);
   fChain->SetBranchAddress("tpc_s1_start_time", &tpc_s1_start_time, &b_tpc_s1_start_time);
   fChain->SetBranchAddress("tpc_total_s1", &tpc_total_s1, &b_tpc_total_s1);
   fChain->SetBranchAddress("tpc_total_s1_f90", &tpc_total_s1_f90, &b_tpc_total_s1_f90);
   fChain->SetBranchAddress("tpc_total_s2", &tpc_total_s2, &b_tpc_total_s2);
   fChain->SetBranchAddress("tpc_total_s2_f90", &tpc_total_s2_f90, &b_tpc_total_s2_f90);
   fChain->SetBranchAddress("tpc_t_drift", &tpc_t_drift, &b_tpc_t_drift);
   fChain->SetBranchAddress("tpc_s1_late", &tpc_s1_late, &b_tpc_s1_late);
   fChain->SetBranchAddress("tpc_npulses", &tpc_npulses, &b_tpc_npulses);
   fChain->SetBranchAddress("tpc_timestamp", &tpc_timestamp, &b_tpc_timestamp);
   fChain->SetBranchAddress("tpc_live_time", &tpc_live_time, &b_tpc_live_time);
   fChain->SetBranchAddress("tpc_inhibit_time", &tpc_inhibit_time, &b_tpc_inhibit_time);
   fChain->SetBranchAddress("tpc_has_s3", &tpc_has_s3, &b_tpc_has_s3);
   fChain->SetBranchAddress("tpc_masa_x", &tpc_masa_x, &b_tpc_masa_x);
   fChain->SetBranchAddress("tpc_masa_y", &tpc_masa_y, &b_tpc_masa_y);
   fChain->SetBranchAddress("tpc_jason_x", &tpc_jason_x, &b_tpc_jason_x);
   fChain->SetBranchAddress("tpc_jaosn_y", &tpc_jaosn_y, &b_tpc_jason_y);
   fChain->SetBranchAddress("od_eventID", &od_eventID, &b_od_eventID);
   fChain->SetBranchAddress("od_gps_fine", &od_gps_fine, &b_od_gps_fine);
   fChain->SetBranchAddress("od_gps_coarse", &od_gps_coarse, &b_od_gps_coarse);
   fChain->SetBranchAddress("od_timestamp", &od_timestamp, &b_od_timestamp);
   fChain->SetBranchAddress("od_nclusters", &od_nclusters, &b_od_nclusters);
   fChain->SetBranchAddress("od_lsv_charge", &od_lsv_charge, &b_od_lsv_charge);
   fChain->SetBranchAddress("od_wt_charge", &od_wt_charge, &b_od_wt_charge);
   fChain->SetBranchAddress("od_cluster_charge", &od_cluster_charge, &b_od_cluster_charge);
   fChain->SetBranchAddress("od_cluster_start", &od_cluster_start, &b_od_cluster_start);
   fChain->SetBranchAddress("od_cluster_height", &od_cluster_height, &b_od_cluster_height);
   fChain->SetBranchAddress("od_cluster_multiplicity", &od_cluster_multiplicity, &b_od_cluster_multiplicity);
   fChain->SetBranchAddress("od_cluster_pass_multcut", &od_cluster_pass_multcut, &b_od_cluster_pass_multcut);
   fChain->SetBranchAddress("od_cluster_dtprompt", &od_cluster_dtprompt, &b_od_cluster_dtprompt);
   fChain->SetBranchAddress("od_roi_lsv_id", &od_roi_lsv_id, &b_od_roi_lsv_id);
   fChain->SetBranchAddress("od_roi_lsv_charge", &od_roi_lsv_charge, &b_od_roi_lsv_charge);
   fChain->SetBranchAddress("od_roi_lsv_max_multiplicity", &od_roi_lsv_max_multiplicity, &b_od_roi_lsv_max_multiplicity);
   fChain->SetBranchAddress("od_slider_lsv_id", &od_slider_lsv_id, &b_od_slider_lsv_id);
   fChain->SetBranchAddress("od_slider_lsv_charge", &od_slider_lsv_charge, &b_od_slider_lsv_charge);
   fChain->SetBranchAddress("od_slider_lsv_time", &od_slider_lsv_time, &b_od_slider_lsv_time);
   fChain->SetBranchAddress("od_slider_lsv_max_multiplicity", &od_slider_lsv_max_multiplicity, &b_od_slider_lsv_max_multiplicity);
}

Bool_t SLADDSTSelector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef SLADDSTSelector_cxx
