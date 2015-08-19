#include "vetopulseslad.hh"

using namespace std;
ClassImp(vetopulseslad);

void vetopulseslad::SLADProcess(int start, int end)
{
  SetSLADTPCIndir(GetRealinputdir());
  SetSLADODIndir(GetRealinputdir());
  cout<<GetSLADTPCIndir()<<"  "<<GetSLADODIndir()<<endl;

  Create_TPCChain();
  Create_ODChain();
  for(int i=start; i<=end; i++){
    if(!Add_TPCTree(i) || !Add_ODTree(i))
      continue;
  }
  Add_TPCFriend();
  tpc_chain->AddFriend(od_chain);
  Load_TPCTree(tpc_chain,e);
  Load_ODTree(tpc_chain,t);    
  BookHistograms();
  LoopOverEvent();
  SaveHistograms();
}

void vetopulseslad::BookHistograms()
{
  int Bins = 2000;
  livetime = new TNtuple("livetime","tpc live time","life_time:eventid");
  FullSpectrum = new TH1F("FullSpectrum","lsv cluster fixed width charge full spectrum;charge [PE];Counts",Bins,0,2000);
  FullSpectrum_ntuple = new TNtuple("FullSpectrum_ntuple","Veto All Data","charge:height:multi:time");
  F1Spectrum = new TH1F("F1Spectrum","TPC energy trigger spectrum;charge [PE];Counts",Bins,0,2000);
  F2Spectrum = new TH1F("F2Spectrum","Late coincidece spectrum;charge [PE];Counts",Bins,0,2000);
  C14Spectrum = new TH1F("C14Spectrum","lsv cluster fixed width charge;charge [PE];Counts",300,0,300);
  Co60Spectrum = new TH1F("Co60Spectrum","Co60 lsv cluster fixed width charge;charge [PE];Counts",1000,0,1000);
  RandomSpectrum = new TH1F("RandomSpectrum","Random spectrum;charge [PE];Counts",Bins,0,2000);
  Veto_TPC_hist = new TH2F("Veto_TPC_hist","lsv cluster charge vs tpc s1 TBAsym;total_s1_TBAcorr [PE];od_cluster_charge[PE]",1000,0,25.e+3,1000,0,2000);
  
}

void vetopulseslad::LoopOverEvent()
{
  int nEntries = tpc_chain->GetEntries();
  cout<<nEntries<<endl;
  Double_t par2[] = { -0.0397956, -0.27216, 0.794036,
		      1.70427, -3.98323, -8.50783, -2.66051 };
  Double_t TBAsym = (e.total_s1_top - e.total_s1_bottom) /
    (e.total_s1_top + e.total_s1_bottom);
  Double_t x = TBAsym; // normalized at 200 PE
  Double_t diff_total_s1 = par2[0] + (par2[1] + (par2[2] +
						 (par2[3] + (par2[4] + par2[5] * x) * x) * x) * x) *x;//(total_s1-total_s1_corr)/total_s1
  Double_t total_s1_TBAcorr = e.s1 * (1. - diff_total_s1);
  
  for(int j=0; j<nEntries; j++){
    tpc_chain->GetEntry(j);
    bool afterpulse = false;
    livetime->Fill(e.life_time,e.event_id);
    if(t.od_present==1){      
      for(size_t i=0; i<t.od_cluster_start->size(); i++){
	if(t.od_cluster_charge->at(i)>0 && t.od_cluster_pass_multcut->at(i)==1 && t.wt_total_charge<500){	
	  double gps = t.od_cluster_dtprompt->at(i);
	  RandomSpectrum->Fill(t.od_cluster_charge->at(i));
	  //  if(TMath::Abs(gps)<0.05)
	  if(gps>-2 && gps<0.1){
	    FullSpectrum->Fill(t.od_cluster_charge->at(i));	
	    FullSpectrum_ntuple->Fill(t.od_cluster_charge->at(i),t.od_cluster_height->at(i),t.od_cluster_multiplicity->at(i),gps);
	    Veto_TPC_hist->Fill(total_s1_TBAcorr,t.od_cluster_charge->at(i));
	    if(total_s1_TBAcorr>10200 && total_s1_TBAcorr<11600 && t.od_cluster_charge->at(i)>400 && t.od_cluster_charge->at(i)<800)
	      Co60Spectrum->Fill(t.od_cluster_charge->at(i));
	  }
	  if((gps<-6 && gps>-10.4)|| gps>30.)	 
	    C14Spectrum->Fill(t.od_cluster_charge->at(i));     
	  if(gps>0.1 && gps<2 && t.od_cluster_charge->at(i)>50)
	    afterpulse = true;	
	}
      }
      for(size_t i=0; !afterpulse && i<t.od_cluster_start->size(); i++){
	if(t.od_cluster_charge->at(i)>0 && t.od_cluster_pass_multcut->at(i)==1 && t.wt_total_charge<500){
	  double gps = t.od_cluster_dtprompt->at(i);
	  F1Spectrum->Fill(t.od_cluster_charge->at(i));
	  if(gps>2.)
	    F2Spectrum->Fill(t.od_cluster_charge->at(i));	
	}
      }
    }  
  }
}

void vetopulseslad::SaveHistograms()
{
  string output = GetRealoutputdir()+GetRealFile();  
  TFile outfile(output.c_str(), "RECREATE");
  livetime->Write();
  RandomSpectrum->Write();
  FullSpectrum->Write();
  FullSpectrum_ntuple->Write();
  Veto_TPC_hist->Write();
  Co60Spectrum->Write();
  C14Spectrum->Write();
  F1Spectrum->Write();
  F2Spectrum->Write();

  outfile.Write();
  outfile.Close();
  cout<<"Successfully Saved Processed SLAD Data."<<endl;
}


/*
bool vetopulseslad::Load_TPCTree(TChain* tpctree, SLADTPCEvent& e)
{
  //events tree
  tpctree->SetBranchAddress("run_id", &e.run_id);
  tpctree->SetBranchAddress("subrun_id", &e.subrun_id);
  tpctree->SetBranchAddress("event_id", &e.event_id);
  tpctree->SetBranchAddress("trigger_multiplicity", &e.trigger_multiplicity);
  tpctree->SetBranchAddress("tpc_digital_sum", &e.tpc_digital_sum);
  tpctree->SetBranchAddress("hasV1724", &e.hasV1724);
  //logbook tree
  tpctree->SetBranchAddress("logbook.integrated_lifetime", &e.integrated_lifetime);
  tpctree->SetBranchAddress("logbook.integrated_inhibittime", &e.integrated_inhibitime);
  //gps tree
  tpctree->SetBranchAddress("gps.gps_coarse", &e.gps_coarse);
  tpctree->SetBranchAddress("gps.gps_fine", &e.gps_fine);
  //nchannel tree
  tpctree->SetBranchAddress("nchannel.nchannel", &e.nchannel);
  //baseline tree
  tpctree->SetBranchAddress("baseline.SumChannelHasNoBaseline", &e.baseline_not_found);
  //long_lifetime tree
  tpctree->SetBranchAddress("long_lifetime.lifetime", &e.life_time);
  tpctree->SetBranchAddress("long_lifetime.inhibittime", &e.inhibit_time);
  //acqui_window tree
  tpctree->SetBranchAddress("acqui_window.acqui_window", &e.acqui_window);
  //standard_cuts tree
  tpctree->SetBranchAddress("standard_cuts.selection_results", &e.selection_results);
  //tdrift tree
  tpctree->SetBranchAddress("tdrift.tdrift", &e.tdrift);
  //npulses tree
  tpctree->SetBranchAddress("npulses.n_phys_pulses", &e.npulses);
  tpctree->SetBranchAddress("npulses.has_s3", &e.has_s3);
  tpctree->SetBranchAddress("npulses.has_s1echo", &e.has_s1echo);
  //s1_time tree
  tpctree->SetBranchAddress("s1_time.s1_start_time", &e.s1_start_time);
  tpctree->SetBranchAddress("s1_time.s1_end_time", &e.s1_end_time);
  tpctree->SetBranchAddress("s1_time.s2_start_time", &e.s2_start_time);
  tpctree->SetBranchAddress("s1_time.s2_end_time", &e.s2_end_time);
  //s1 tree
  tpctree->SetBranchAddress("s1.total_s1", &e.s1);
  tpctree->SetBranchAddress("s1.total_s1_corr", &e.s1_corr);
  tpctree->SetBranchAddress("s1.total_s1_top", &e.total_s1_top);
  tpctree->SetBranchAddress("s1.total_s1_bottom", &e.total_s1_bottom);
  //s1_f90 tree
  tpctree->SetBranchAddress("s1_f90.total_f90", &e.s1_total_f90);
  tpctree->SetBranchAddress("s1_f90.total_f90_fixed", &e.s1_total_f90_fixed);
  tpctree->SetBranchAddress("s1_f90.total_f90_spe_mean", &e.s1_total_f90_spe_mean);
  //s1_saturation tree
  tpctree->SetBranchAddress("s1_saturation.is_saturated_pulse0", &e.is_saturated_pulse0);
  //s1_fraction tree
  tpctree->SetBranchAddress("s1_fraction.s1_max_chan", &e.s1_max_chan);
  tpctree->SetBranchAddress("s1_fraction.s1_max_frac", &e.s1_max_frac);
  //max_s1_frac_cut tree
  tpctree->SetBranchAddress("max_s1_frac_cut.max_s1_frac_cut_threshold99", &e.max_s1_frac_cut_threshold99);
  tpctree->SetBranchAddress("max_s1_frac_cut.max_s1_frac_cut_exceeds99", &e.max_s1_frac_cut_exceeds99);
  //s2 tree
  tpctree->SetBranchAddress("s2.total_s2", &e.s2);
  tpctree->SetBranchAddress("s2.total_s2_corr", &e.s2_corr);
  tpctree->SetBranchAddress("s2.total_s2_top", &e.total_s2_top);
  tpctree->SetBranchAddress("s2.total_s2_bottom", &e.total_s2_bottom);
  //s2_f90 tree
  tpctree->SetBranchAddress("s2_f90.total_s2_f90", &e.s2_total_f90);
  tpctree->SetBranchAddress("s2_f90.total_s2_f90_fixed", &e.s2_total_f90_fixed);
  tpctree->SetBranchAddress("s2_f90.total_s2_f90_spe_mean", &e.s2_total_f90_spe_mean);
  //s2_saturation tree
  tpctree->SetBranchAddress("s2_saturation.is_saturated_pulse1", &e.is_saturated_pulse1);
  //bary_s2 tree
  tpctree->SetBranchAddress("bary_s2.bary_s2_x", &e.bary_s2_x);
  tpctree->SetBranchAddress("bary_s2.bary_s2_y", &e.bary_s2_y);
  //s2_fraction tree
  tpctree->SetBranchAddress("s2_fraction.s2_chan", e.s2_chan);
  tpctree->SetBranchAddress("s2_fraction.s2_max_chan", &e.s2_max_chan);
  tpctree->SetBranchAddress("s2_fraction.s2_max_frac", &e.s2_max_frac);
  //xylocator_xy tree
  tpctree->SetBranchAddress("xylocator_xy.xyl_SCM", &e.jason_SCM);
  tpctree->SetBranchAddress("xylocator_xy.xyl_best_x", &e.jason_x);
  tpctree->SetBranchAddress("xylocator_xy.xyl_best_y", &e.jason_y);
  tpctree->SetBranchAddress("xylocator_xy.xyl_best_chi2", &e.jason_chi2);
  tpctree->SetBranchAddress("xylocator_xy.xyl_best_r", &e.jason_r);
  tpctree->SetBranchAddress("xylocator_xy.xyl_best_theta", &e.jason_theta);
  tpctree->SetBranchAddress("xylocator_xy.xyl_best_xycorr_factor", &e.jason_xycorr_factor);
  //allpulses_xyl_xy tree
  tpctree->SetBranchAddress("allpulses_xyl_xy.allpulses_xyl_npulses", &e.allpulses_xyl_npulses);
  tpctree->SetBranchAddress("allpulses_xyl_xy.allpulses_xyl_fidcut_ratio", &e.allpulses_xyl_fidcut_ratio);
  tpctree->SetBranchAddress("allpulses_xyl_xy.allpulses_xyl_x", &e.allpulses_xyl_x);
  tpctree->SetBranchAddress("allpulses_xyl_xy.allpulses_xyl_y", &e.allpulses_xyl_y);
  tpctree->SetBranchAddress("allpulses_xyl_xy.allpulses_xyl_chi2", &e.allpulses_xyl_chi2);
  tpctree->SetBranchAddress("allpulses_xyl_xy.allpulses_xyl_r", &e.allpulses_xyl_r);
  tpctree->SetBranchAddress("allpulses_xyl_xy.allpulses_xyl_theta", &e.allpulses_xyl_theta);
  tpctree->SetBranchAddress("allpulses_xyl_xy.allpulses_xyl_xycorr_factor", &e.allpulses_xyl_xycorr_factor);
  //masa_xy tree
  tpctree->SetBranchAddress("masas_xy.masas_x", &e.masa_x);
  tpctree->SetBranchAddress("masas_xy.masas_y", &e.masa_y);
  tpctree->SetBranchAddress("masas_xy.masas_chi2", &e.masa_chi2);
  tpctree->SetBranchAddress("masas_xy.masas_r", &e.masa_r);
  tpctree->SetBranchAddress("masas_xy.masas_theta", &e.masa_theta);
  tpctree->SetBranchAddress("masas_xy.masas_xycorr_factor", &e.masa_xycorr_factor);
  //allpulses_xy tree
  tpctree->SetBranchAddress("allpulses_xy.allpulses_npulses", &e.allpulses_npulses);
  tpctree->SetBranchAddress("allpulses_xy.allpulses_x", &e.allpulses_x);
  tpctree->SetBranchAddress("allpulses_xy.allpulses_y", &e.allpulses_y);
  tpctree->SetBranchAddress("allpulses_xy.allpulses_chi2", &e.allpulses_chi2);
  tpctree->SetBranchAddress("allpulses_xy.allpulses_r", &e.allpulses_r);
  tpctree->SetBranchAddress("allpulses_xy.allpulses_theta", &e.allpulses_theta);
  tpctree->SetBranchAddress("allpulses_xy.allpulses_xycorr_factor", &e.allpulses_xycorr_factor);
  //pulse_info tree
  tpctree->SetBranchAddress("pulse_info.pulse_info_npulses",&e.pulse_info_npulses);
  tpctree->SetBranchAddress("pulse_info.pulse_info_pulse_id",&e.pulse_info_pulse_id);
  tpctree->SetBranchAddress("pulse_info.pulse_info_start_time",&e.pulse_info_start_time);
  tpctree->SetBranchAddress("pulse_info.pulse_info_end_time",&e.pulse_info_end_time);
  tpctree->SetBranchAddress("pulse_info.pulse_info_peak_time",&e.pulse_info_peak_time);
  tpctree->SetBranchAddress("pulse_info.pulse_info_max_chan",&e.pulse_info_max_chan);
  tpctree->SetBranchAddress("pulse_info.pulse_info_total_npe",&e.pulse_info_total_npe);
  tpctree->SetBranchAddress("pulse_info.pulse_info_f90",&e.pulse_info_f90);
  tpctree->SetBranchAddress("pulse_info.pulse_info_f200",&e.pulse_info_f200);
  tpctree->SetBranchAddress("pulse_info.pulse_info_saturated",&e.pulse_info_saturated);
  tpctree->SetBranchAddress("pulse_info.pulse_info_fixed_int1",&e.pulse_info_fixed_int1);
  tpctree->SetBranchAddress("pulse_info.pulse_info_fixed_int2",&e.pulse_info_fixed_int2);
  tpctree->SetBranchAddress("pulse_info.pulse_info_ch_light_size",&e.pulse_info_ch_light_size);
  tpctree->SetBranchAddress("pulse_info.pulse_info_ch_light",&e.pulse_info_light);
  tpctree->SetBranchAddress("pulse_info.pulse_info_satcorr_f90",&e.pulse_info_satcorr_f90);
  tpctree->SetBranchAddress("pulse_info.pulse_info_satcorr_f90_fixed",&e.pulse_info_satcorr_f90_fixed);
  tpctree->SetBranchAddress("pulse_info.pulse_info_satcorr_fixed_int1",&e.pulse_info_satcorr_fixed_int1);
  tpctree->SetBranchAddress("pulse_info.pulse_info_satcorr_fixed_int2",&e.pulse_info_satcorr_fixed_int2);
  tpctree->SetBranchAddress("pulse_info.pulse_info_satcorr_f200",&e.pulse_info_satcorr_f200);
  tpctree->SetBranchAddress("pulse_info.pulse_info_satcorr_f200_fixed",&e.pulse_info_satcorr_f200_fixed);

  return true;
}

bool vetopulseslad::Load_ODTree(TChain* odtree, SLADODEvent &t)
{
  odtree->SetBranchAddress("veto_run_id",&t.od_run);
  odtree->SetBranchAddress("veto_event_id",&t.od_event);
  odtree->SetBranchAddress("veto_present",&t.od_present);
  odtree->SetBranchAddress("veto_timestamp_us", &t.od_timestamp);
  odtree->SetBranchAddress("veto_lsv_total_charge", &t.lsv_total_charge);
  odtree->SetBranchAddress("veto_wt_total_charge", &t.wt_total_charge);
  odtree->SetBranchAddress("veto_roi_lsv_charge_vec",&t.roi_lsv_charge);
  odtree->SetBranchAddress("veto_slider_lsv_charge_vec",&t.slider_lsv_charge);
  odtree->SetBranchAddress("veto_slider_lsv_time_vec",&t.slider_lsv_time);

  odtree->SetBranchAddress("veto_nclusters", &t.od_nclusters);
  odtree->SetBranchAddress("veto_cluster_charge", &t.od_cluster_charge);
  odtree->SetBranchAddress("veto_cluster_start", &t.od_cluster_start);
  odtree->SetBranchAddress("veto_cluster_height", &t.od_cluster_height);
  odtree->SetBranchAddress("veto_cluster_multiplicity", &t.od_cluster_multiplicity);
  odtree->SetBranchAddress("veto_cluster_pass_multcut", &t.od_cluster_pass_multcut);
  odtree->SetBranchAddress("veto_cluster_dtprompt", &t.od_cluster_dtprompt);

  return true;
}

void vetopulseslad::Create_TPCChain()
{
  tpc_chain = new TChain("events");
  s2_fraction = new TChain("s2_fraction");
  masaxy = new TChain("masas_xy");
  jasonxy = new TChain("xylocator_xy");
  pulse_info = new TChain("pulse_info");
}

bool vetopulseslad::Add_TPCTree(int i)
{
  TString tpcfile;
  tpcfile.Form("Run%06d.root",i);
  tpcfile.Prepend(GetRealinputdir());

  if(!Verifydatafile(tpcfile))
    return false;

  TString s2file = tpcfile;
  s2file.Remove(s2file.Length()-5);
  s2file+="_s2.root";

  TString masaxyfile = tpcfile;
  masaxyfile.Remove(masaxyfile.Length()-5);
  masaxyfile+="_masas_xy.root";

  TString jasonxyfile = tpcfile;
  jasonxyfile.Remove(jasonxyfile.Length()-5);
  jasonxyfile+="_xylocator_xy.root";

  TString pulsefile = tpcfile;
  pulsefile.Remove(pulsefile.Length()-5);
  pulsefile+="_allpulses.root";

  if(!Verifydatafile(s2file) || !Verifydatafile(masaxyfile) || !Verifydatafile(jasonxyfile) || !Verifydatafile(pulsefile))
    return false;

  tpc_chain->Add(tpcfile);
  s2_fraction->Add(s2file);
  masaxy->Add(masaxyfile);
  jasonxy->Add(jasonxyfile);
  pulse_info->Add(pulsefile);
  return true;
}

void vetopulseslad::Add_TPCFriend()
{
  tpc_chain->AddFriend(s2_fraction);
  tpc_chain->AddFriend(masaxy);
  tpc_chain->AddFriend(jasonxy);
  tpc_chain->AddFriend(pulse_info);
}

void vetopulseslad::Create_ODChain()
{
  od_chain = new TChain("veto");
}

bool vetopulseslad::Add_ODTree(int i)
{
  TString odfile;
  odfile.Form("Run%06d_veto.root",i);
  odfile.Prepend(GetRealinputdir());
  if(!Verifydatafile(odfile))
    return false;
  else
    od_chain->Add(odfile);
  return true;
}
*/
