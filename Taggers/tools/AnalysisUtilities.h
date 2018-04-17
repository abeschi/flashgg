#ifndef ANALYSYS_UTILITIES_H
#define ANALYSYS_UTILITIES_H




#ifndef CMS_LUMI_H
#include "CMS_lumi.h"
#endif

#ifndef CMS_STYLE
#include "tdrstyle.h"
#endif

#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPaveStats.h"
#include "TLegend.h"
#include "TChain.h"
#include "TVirtualFitter.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TMath.h"
#include "TBox.h"
#include "TLorentzVector.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;




float compute_R(float eta_1, float eta_2, float phi_1, float phi_2);

void MakePlot(TH1F* histo, TString title);

std::vector<float>  parse (TString fileName);

bool bjetCut(int loose, int medium, int tight, int loose_cut, int medium_cut, int tight_cut);

std::vector<float> FindSmallestInterval(TH1F* histo, const float& fraction=0.68, const bool& verbosity=0);

TF1* FitBkg(TH1F* histo);

void packageResults(TH1F* Signal, TH1F* Bkg, TH1F* Cs, TF1* fitFunc, float min, float max, std::vector<float> thresholds, float* purity, float lumiFactor, TString title);










#endif




     






