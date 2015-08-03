#ifndef vetoreadcfg_H
#define vetoreadcfg_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

struct Config{
  string mcindir;
  string mcfile;
  string realindir;
  string realfile;
  string outdir;
  string outfile;
  int startrun;
  int endrun;
  int mc_tag;
  int real_tag;
  int totalfit_tag;
  int startfit;
  bool DST;
  bool OD;
  bool SLADDST;
  bool fieldon;
  
  void loadConfig(char*);

};

//void loadConfig(Config& config,char* cfile) {
inline void Config::loadConfig(char* cfile) {
  ifstream fin(cfile);
  string line;
  while (getline(fin, line)) {
    istringstream sin(line.substr(line.find("=") + 1));
    if ((int)line.find("mcindir") != -1)
      sin >> mcindir;
    else if ((int)line.find("mcfile") != -1)
      sin >> mcfile;
    else if ((int)line.find("realindir") != -1)
      sin >> realindir;
    else if ((int)line.find("realfile") != -1)
      sin >> realfile;
    else if ((int)line.find("outdir") != -1)
      sin >> outdir;
    else if ((int)line.find("outfile") != -1)
      sin >> outfile;
    else if ((int)line.find("startrun") != -1)
      sin >> startrun;    
    else if ((int)line.find("endrun") != -1)
      sin >> endrun;    
    else if ((int)line.find("mc_tag") != -1)
      sin >> mc_tag;    
    else if ((int)line.find("real_tag") != -1)
      sin >> real_tag;    
    else if ((int)line.find("totalfit_tag") != -1)
      sin >> totalfit_tag;    
    else if ((int)line.find("startfit") != -1)
      sin >> startfit;    
    else if ((int)line.find("DST") != -1)
      sin >> DST;    
    else if ((int)line.find("OD") != -1)
      sin >> OD;    
    else if ((int)line.find("SLADDST") != -1)
      sin >> SLADDST;    
    else if ((int)line.find("fieldon") != -1)
      sin >> fieldon;    

    sin.clear();
  }
}
#endif
