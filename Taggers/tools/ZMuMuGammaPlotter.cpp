/* To compile
g++ -Wall -o ZMuMuGammaPlotter `root-config --cflags --glibs` -L $ROOTSYS/lib -lFoam -lMinuit -lMathMore CMS_lumi.C tdrstyle.C ZMuMuGammaPlotter.cpp
*/


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
#include "TObjArray.h"
#include "TMinuit.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "RoccoR.cc"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <time.h>   


using namespace std;

float compute_DeltaPhi(float phi_1, float phi_2)
{
	float delta_phi = fabs(phi_1 - phi_2);
	if(delta_phi > 3.14159265359)
		delta_phi = 2*3.14159265359 - delta_phi;
	if(phi_1 - phi_2 >= 0)	
		return delta_phi;
	else
		return -delta_phi;
}



void HistoPlotter(TH1F* MC, TH1F* Data, TString name)
{
	TCanvas* c = new TCanvas();

	Data -> SetMarkerStyle(20);
	Data -> SetMarkerSize(1);
	MC -> SetFillStyle(1001);
	MC->SetFillColor(kYellow);

	MC -> Scale((Data->GetSumOfWeights()/MC->GetSumOfWeights()));

	double max = Data->GetBinContent(Data->GetMaximumBin());
	if(MC->GetBinContent(MC->GetMaximumBin()) > max )
		max = MC->GetBinContent(MC->GetMaximumBin());

	MC-> GetYaxis() -> SetRangeUser(0, 1.2*max);
	MC-> GetYaxis() -> SetTitleSize(0.05);
	MC-> GetYaxis() -> SetTitleFont(42);
	MC-> GetYaxis() -> SetLabelSize(0.045);
	MC-> GetYaxis() -> SetLabelFont(42);
	MC-> GetXaxis() -> SetLabelSize(0);

	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 0.97);
	pad1->SetBottomMargin(0.035);
	pad1->Draw();
	pad1->cd();

	MC -> Draw("histo");
	Data -> Draw("EP SAME");


	TLegend* leg = new TLegend(0.2, 0.80, 0.29, 0.92);
	leg -> AddEntry(Data, "Data", "p");
	leg -> AddEntry(MC, "MC", "f");
	leg -> Draw("SAME");


	c->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.10, 1, 0.35);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0);
	pad2 -> SetGridy();
	pad2->Draw();
	pad2->cd();

	TH1F *h = (TH1F*)Data->Clone("h");
	h->SetLineColor(kBlack);

	if(name=="Mass2" || name=="MassEB2" || name=="MassEE2")
	{
		h->SetMinimum(0.);  // Define Y ..
		h->SetMaximum(2.); // .. range
	}
	else
	{
		h->SetMinimum(0.5);  // Define Y ..
		h->SetMaximum(1.5); // .. range
	}

	h->Sumw2();
	h->SetStats(0);      // No statistics on lower plot
	h->Divide(MC);
	h->SetMarkerStyle(21);
	h -> SetTitle("");
	h-> GetYaxis() -> SetTitle("Data/MC");
	//h-> GetYaxis() -> SetRangeUser(0., 2.);


	// Y axis ratio plot settings
	h->GetYaxis()->SetNdivisions(-10);
	h->GetYaxis()->SetTitleSize(0.13);
	h->GetYaxis()->SetTitleFont(42);
	h->GetYaxis()->SetTitleOffset(0.5);
	h->GetYaxis()->SetLabelFont(42); // Absolute font size in pixel (precision 3)
	h->GetYaxis()->SetLabelSize(0.12);

	// X axis ratio plot settings
	h->GetXaxis()->SetTitleSize(0.15);
	h->GetXaxis()->SetTitleFont(42);
	h->GetXaxis()->SetLabelFont(42); // Absolute font size in pixel (precision 3)
	h->GetXaxis()->SetLabelSize(0.12);

	h->Draw("EP");

	CMS_lumi(c, 0, 0);


	c -> SaveAs(name + ".png");
	c -> SaveAs(name + ".pdf");

	delete c;
	delete leg;

	return;
}




void HistoPlotter2D(TH2F* h, TString name)
{
	TCanvas* c = new TCanvas();

	h -> GetZaxis() -> SetTitleSize(0.05);
	h -> GetZaxis() -> SetTitleOffset(1.5);
	h -> GetZaxis() -> SetLabelSize(0.04);
	h -> Draw("COLZ");
	CMS_lumi(c, 0, 0);

	c -> SetRightMargin(0.2);
	c -> SetBottomMargin(0.15);

	c -> SaveAs(name + ".png");
	c -> SaveAs(name + ".pdf");

	delete c;
	return;
}


int main(int argc, char *argv[])
{
	writeExtraText = true;       // if extra text
	extraText  = "Preliminary";  // default extra text is "Preliminary"
	lumi_sqrtS = "36.8 fb^{-1} (13 TeV)";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

	setTDRStyle();
	gStyle -> SetOptFit(0);
	gStyle -> SetOptStat(0);

	srand (time(NULL));

	bool DoNotApplyScaleAndSmearing = 0;
	//Final plots

	TH1F* MCMass = new TH1F("MCMass", "; M_{Z} (GeV); Counts", 40 , 80, 100);
	TH1F* MCMassEB = new TH1F("MCMassEB", "; M_{Z} (GeV); Counts", 40 , 80, 100);
	TH1F* MCMassEE = new TH1F("MCMassEE", "; M_{Z} (GeV); Counts", 40 , 80, 100);
	TH1F* MCMass2 = new TH1F("MCMass2", "; M_{Z} (GeV); Counts", 60 , 60, 120);
	TH1F* MCMassEB2 = new TH1F("MCMassEB2", "; M_{Z} (GeV); Counts", 60 , 60, 120);
	TH1F* MCMassEE2 = new TH1F("MCMassEE2", "; M_{Z} (GeV); Counts", 60 , 60, 120);
	TH1F* MCMassEBHighR9 = new TH1F("MCMassEBHighR9", "; M_{Z} (GeV); Counts", 40 , 80, 100);
	TH1F* MCMassEEHighR9 = new TH1F("MCMassEEHighR9", "; M_{Z} (GeV); Counts", 40 , 80, 100);
	TH1F* MCMassEBLowR9 = new TH1F("MCMassEBLowR9", "; M_{Z} (GeV); Counts", 40 , 80, 100);
	TH1F* MCMassEELowR9 = new TH1F("MCMassEELowR9", "; M_{Z} (GeV); Counts", 40 , 80, 100);
	TH1F* MCMuMuMass = new TH1F("MCMuMuMass", "; M_{#mu#mu} (GeV); Counts", 70 , 30, 100);
	TH1F* MCMuMuPt = new TH1F("MCMuMuPt", "; P_{T #mu#mu} (GeV); Counts", 50 , 0, 100);
	TH1F* MCLeadingMuPt = new TH1F("MCLeadingMuPt", "; P_{T} (GeV); Counts",50 , 0, 100);
	TH1F* MCLeadingMuEta = new TH1F("MCLeadingMuEta", "; #eta; Counts", 50 , -3, 3);
	TH1F* MCLeadingMuPhi = new TH1F("MCLeadingMuPhi", "; #varphi; Counts", 50 , -3.15, 3.15);
	TH1F* MCLeadingMuCharge = new TH1F("MCLeadingCharge", "; Charge; Counts", 10 , -2, 2);
	TH1F* MCLeadingMuNtrk = new TH1F("MCLeadingMuNtrk", "; N_{tkr}; Counts", 31 , -0.5, 30.5);        // For Rochester corrections
	TH1F* MCSubleadingMuPt = new TH1F("MCSubleadingMuPt", "; P_{T} (GeV); Counts", 50 , 0, 50);
	TH1F* MCSubleadingMuEta = new TH1F("MCSubleadingMuEta", "; #eta; Counts", 50 , -3, 3);
	TH1F* MCSubleadingMuPhi = new TH1F("MCSubleadingMuPhi", "; #varphi; Counts", 50 , -3.15, 3.15);
	TH1F* MCSunleadingMuCharge = new TH1F("MCSubleadingCharge", "; Charge; Counts", 10 , -2, 2);
	TH1F* MCSubleadingMuNtrk = new TH1F("MCSubleadingMuNtrk", "; N_{tkr}; Counts",  31 , -0.5, 30.5);        // For Rochester corrections
	TH1F* MCPhotonE = new TH1F("MCPhotonE", "; Energy (GeV); Counts", 45 , 10, 100);
	TH1F* MCPhotonPt = new TH1F("MCPhotonPt", "; P_{T} (GeV); Counts", 60, 10, 70);
	TH1F* MCPhotonEta = new TH1F("MCPhotonEta", "; #eta; Counts", 50 , -2.5, 2.5);
	TH1F* MCPhotonPhi = new TH1F("MCPhotonPhi", "; #varphi; Counts", 50 , -3.15, 3.15);
	TH1F* MCPhotonR9 = new TH1F("MCPhotonR9", "; R_{9}; Counts", 50 , 0.5, 1);
	TH1F* MCPhotonR9EB = new TH1F("MCPhotonR9EB", "; R_{9}; Counts", 50 , 0.5, 1);
	TH1F* MCPhotonR9EE = new TH1F("MCPhotonR9EE", "; R_{9}; Counts", 50 , 0.5, 1);
	TH1F* MCPhotonS4 = new TH1F("MCPhotonS4", "; S_{4}; Counts", 50 , 0.5, 1);
	TH1F* MCPhotonS4EB = new TH1F("MCPhotonS4EB", "; S_{4}; Counts", 50 , 0.5, 1);
	TH1F* MCPhotonS4EE = new TH1F("MCPhotonS4EE", "; S_{4}; Counts", 50 , 0.5, 1);
	TH1F* MCPhotonEtaWidth = new TH1F("MCPhotonEtaWidth", "; EtaWidth; Counts", 50 , 0., 0.04);
	TH1F* MCPhotonEtaWidthEB = new TH1F("MCPhotonEtaWidthEB", "; EtaWidth; Counts", 50 , 0., 0.02);
	TH1F* MCPhotonEtaWidthEE = new TH1F("MCPhotonEtaWidthEE", "; EtaWidth; Counts", 50 , 0., 0.04);
	TH1F* MCPhotonSigmaIetaIeta = new TH1F("MCPhotonSigmaIEtaIEta", "; #sigma_{I#etaI#eta}; Counts", 50, 0.005, 0.035);
	TH1F* MCPhotonSigmaIetaIetaEB = new TH1F("MCPhotonSigmaIEtaIEtaEB", "; #sigma_{I#etaI#eta}; Counts", 50, 0.006, 0.012);
	TH1F* MCPhotonSigmaIetaIetaEE = new TH1F("MCPhotonSigmaIEtaIEtaEE", "; #sigma_{I#etaI#eta}; Counts", 50, 0.020, 0.035);
	TH1F* MCPhotonSigmaEtaEta = new TH1F("MCPhotonSigmaEtaEta", "; #sigma_{#eta#eta}; Counts", 70 , 0.005, 0.035);
	TH1F* MCPhotonSigmaEtaEtaEB = new TH1F("MCPhotonSigmaEtaEtaEB", "; #sigma_{#eta#eta}; Counts", 70 , 0.006, 0.012);
	TH1F* MCPhotonSigmaEtaEtaEE = new TH1F("MCPhotonSigmaEtaEtaEE", "; #sigma_{#eta#eta}; Counts", 70 , 0.005, 0.035);

	TH1F* DataMass = new TH1F("DataMass", "; M_{Z} (GeV); Counts", 40 , 80, 100);
	TH1F* DataMassEB = new TH1F("DataMassEB", "; M_{Z} (GeV); Counts", 40 , 80, 100);
	TH1F* DataMassEE = new TH1F("DataMassEE", "; M_{Z} (GeV); Counts", 40 , 80, 100);
	TH1F* DataMass2 = new TH1F("DataMass2", "; M_{Z} (GeV); Counts", 60 , 60, 120);
	TH1F* DataMassEB2 = new TH1F("DataMassEB2", "; M_{Z} (GeV); Counts", 60 , 60, 120);
	TH1F* DataMassEE2 = new TH1F("DataMassEE2", "; M_{Z} (GeV); Counts", 60 , 60, 120);
	TH1F* DataMassEBHighR9 = new TH1F("DataMassEBHighR9", "; M_{Z} (GeV); Counts", 40 , 80, 100);
	TH1F* DataMassEEHighR9 = new TH1F("DataMassEEHighR9", "; M_{Z} (GeV); Counts", 40 , 80, 100);
	TH1F* DataMassEBLowR9 = new TH1F("DataMassEBLowR9", "; M_{Z} (GeV); Counts", 40 , 80, 100);
	TH1F* DataMassEELowR9 = new TH1F("DataMassEELowR9", "; M_{Z} (GeV); Counts", 40 , 80, 100);
	TH1F* DataMuMuMass = new TH1F("DataMuMuMass", "; M_{#mu#mu} (GeV); Counts", 70 , 30, 100);
	TH1F* DataMuMuPt = new TH1F("DataMuMuPt", "; P_{T #mu#mu} (GeV); Counts", 50 , 0, 100);
	TH1F* DataLeadingMuPt = new TH1F("DataLeadingMuPt", "; P_{T} (GeV); Counts", 50 , 0, 100);
	TH1F* DataLeadingMuEta = new TH1F("DataLeadingMuEta", "; #eta; Counts", 50 , -3, 3);
	TH1F* DataLeadingMuPhi = new TH1F("DataLeadingMuPhi", "; #varphi; Counts", 50 , -3.15, 3.15);
	TH1F* DataLeadingMuCharge = new TH1F("DataLeadingCharge", "; Charge; Counts", 10 , -2, 2);
	TH1F* DataLeadingMuNtrk = new TH1F("DataLeadingMuNtrk", "; N_{tkr}; Counts", 31 , -0.5, 30.5);        // For Rochester corrections
	TH1F* DataSubleadingMuPt = new TH1F("DataSubleadingMuPt", "; P_{T} (GeV); Counts", 50 , 0, 50);
	TH1F* DataSubleadingMuEta = new TH1F("DataSubleadingMuEta", "; #eta; Counts", 50 , -3, 3);
	TH1F* DataSubleadingMuPhi = new TH1F("DataSubleadingMuPhi", "; #varphi; Counts", 50 , -3.15, 3.15);
	TH1F* DataSunleadingMuCharge = new TH1F("DataSubleadingCharge", "; Charge; Counts", 10 , -2, 2);
	TH1F* DataSubleadingMuNtrk = new TH1F("DataSubleadingMuNtrk", "; N_{tkr}; Counts", 31 , -0.5, 30.5);        // For Rochester corrections
	TH1F* DataPhotonE = new TH1F("DataPhotonE", "; Energy (GeV); Counts", 45 , 10, 100);
	TH1F* DataPhotonPt = new TH1F("DataPhotonPt", "; P_{T} (GeV); Counts", 60, 10, 70);
	TH1F* DataPhotonEta = new TH1F("DataPhotonEta", "; #eta; Counts", 50 , -2.5, 2.5);
	TH1F* DataPhotonPhi = new TH1F("DataPhotonPhi", "; #varphi; Counts", 50 , -3.15, 3.15);
	TH1F* DataPhotonR9 = new TH1F("DataPhotonR9", "; R_{9}; Counts", 50 , 0.5, 1);
	TH1F* DataPhotonR9EB = new TH1F("DataPhotonR9EB", "; R_{9}; Counts", 50 , 0.5, 1);
	TH1F* DataPhotonR9EE = new TH1F("DataPhotonR9EE", "; R_{9}; Counts", 50 , 0.5, 1);
	TH1F* DataPhotonS4 = new TH1F("DataPhotonS4", "; S_{4}; Counts", 50, 0.5, 1);
	TH1F* DataPhotonS4EB = new TH1F("DataPhotonS4EB", "; S_{4}; Counts", 50, 0.5, 1);
	TH1F* DataPhotonS4EE = new TH1F("DataPhotonS4EE", "; S_{4}; Counts", 50, 0.5, 1);
	TH1F* DataPhotonEtaWidth = new TH1F("DataPhotonEtaWidth", "; EtaWidth; Counts", 50 , 0., 0.04);
	TH1F* DataPhotonEtaWidthEB = new TH1F("DataPhotonEtaWidthEB", "; EtaWidth; Counts", 50 , 0., 0.02);
	TH1F* DataPhotonEtaWidthEE = new TH1F("DataPhotonEtaWidthEE", "; EtaWidth; Counts", 50 , 0., 0.04);
	TH1F* DataPhotonSigmaIetaIeta = new TH1F("DataPhotonSigmaIEtaIEta", "; #sigma_{I#etaI#eta}; Counts", 50 , 0.005, 0.035);
	TH1F* DataPhotonSigmaIetaIetaEB = new TH1F("DataPhotonSigmaIEtaIEtaEB", "; #sigma_{I#etaI#eta}; Counts", 50 , 0.006, 0.012);
	TH1F* DataPhotonSigmaIetaIetaEE = new TH1F("DataPhotonSigmaIEtaIEtaEE", "; #sigma_{I#etaI#eta}; Counts", 50 , 0.020, 0.035);
	TH1F* DataPhotonSigmaEtaEta = new TH1F("DataPhotonSigmaEtaEta", "; #sigma_{#eta#eta}; Counts", 70 , 0.005, 0.035);
	TH1F* DataPhotonSigmaEtaEtaEB = new TH1F("DataPhotonSigmaEtaEtaEB", "; #sigma_{#eta#eta}; Counts", 70 , 0.006, 0.012);
	TH1F* DataPhotonSigmaEtaEtaEE = new TH1F("DataPhotonSigmaEtaEtaEE", "; #sigma_{#eta#eta}; Counts", 70 , 0.005, 0.035);

	TH2F* MassMuMuVsZ = new TH2F("MassMuMuVsZ", "; M_{#mu#mu} (GeV);  M_{Z} (GeV); Counts", 55, 35, 90, 40, 80, 100);
	TH2F* MassMuMuVsPh = new TH2F("MassMuMuVsPh", "; M_{#mu#mu} (GeV);  p_{T}^{#gamma} (GeV); Counts", 55, 35, 90, 40, 20, 60);
	TH2F* MassTest = new TH2F("MassTest", "; M_{#mu#mu} (GeV);  p_{T}^{#gamma} (GeV); m_{Z}", 55, 35, 90, 40, 20, 60);

	TString MCFolder;
	TString DataFolder;

	if(DoNotApplyScaleAndSmearing)
	{
		MCFolder = "/afs/cern.ch/work/a/abeschi/MuMuGammaMCNoScale/";
		DataFolder = "/afs/cern.ch/work/a/abeschi/MuMuGammaDataNoScale/";
	}

	else
	{
		MCFolder = "/afs/cern.ch/work/a/abeschi/MuMuGammaMC4/";
		DataFolder = "/afs/cern.ch/work/a/abeschi/MuMuGammaData4/";
	}
	int nentries;

	float weight = 0;
	float mass_mmg = 0;
	float mass_mumu = 0;
	float pt_mumu = 0;
	float leadMuonP = 0;
	float leadMuonPt = 0;
	float leadMuonEta = 0;
	float leadMuonPhi = 0;
	float leadMuonCharge = 0;
	float leadMuonNtrk = 0;
	float subleadMuonP = 0;
	float subleadMuonPt = 0;
	float subleadMuonEta = 0;
	float subleadMuonPhi = 0;
	float subleadMuonCharge = 0;
	float subleadMuonNtrk = 0;
	float photonE = 0;
	float photonPt = 0;
	float photonEta = 0;
	float photonPhi = 0;
	float photonR9 = 0;
	float photonS4 = 0;
	float photonEtaWidth = 0;
	float photonSigmaIetaIeta = 0;
	float photonSigmaEtaEta = 0;

//	float muMass = 0.1057;
	int evData = 0;
	int evMC = 0;
	float sumW = 0;

	TChain* DataTree;
	TChain* MCTree;

	float lumiFactor = 36.8;

	DataTree = new TChain("mumugammaDumper/trees/tree");
	DataTree -> Add(DataFolder + "MuMuGammaData*.root");		
	MCTree = new TChain("mumugammaDumper/trees/tree");
	MCTree -> Add(MCFolder + "MuMuGammaMC_*.root");

	MCTree -> SetBranchAddress("weight", &weight);
       	MCTree -> SetBranchAddress("mass_mmg", &mass_mmg);
	MCTree -> SetBranchAddress("mass_mumu", &mass_mumu);
	MCTree -> SetBranchAddress("pt_mumu", &pt_mumu);
	MCTree -> SetBranchAddress("leadMuonP", &leadMuonP);
	MCTree -> SetBranchAddress("leadMuonPt", &leadMuonPt);
	MCTree -> SetBranchAddress("leadMuonEta", &leadMuonEta);
	MCTree -> SetBranchAddress("leadMuonPhi", &leadMuonPhi);
	MCTree -> SetBranchAddress("leadMuonCharge", &leadMuonCharge);
	MCTree -> SetBranchAddress("leadMuonNtrk", &leadMuonNtrk);
	MCTree -> SetBranchAddress("subleadMuonP", &subleadMuonP);
	MCTree -> SetBranchAddress("subleadMuonPt", &subleadMuonPt);
	MCTree -> SetBranchAddress("subleadMuonEta", &subleadMuonEta);
	MCTree -> SetBranchAddress("subleadMuonPhi", &subleadMuonPhi);
	MCTree -> SetBranchAddress("subleadMuonCharge", &subleadMuonCharge);
	MCTree -> SetBranchAddress("subleadMuonNtrk", &subleadMuonNtrk);
	MCTree -> SetBranchAddress("photonE", &photonE);
	MCTree -> SetBranchAddress("photonPt", &photonPt);
	MCTree -> SetBranchAddress("photonEta", &photonEta);
	MCTree -> SetBranchAddress("photonPhi", &photonPhi);
	MCTree -> SetBranchAddress("photonR9", &photonR9);
	MCTree -> SetBranchAddress("photonS4", &photonS4);
	MCTree -> SetBranchAddress("photonEtaWidth", &photonEtaWidth);
	MCTree -> SetBranchAddress("photonSigmaIEtaIEta", &photonSigmaIetaIeta);
	MCTree -> SetBranchAddress("photonEtaEta", &photonSigmaEtaEta);

	DataTree -> SetBranchAddress("mass_mmg", &mass_mmg);
	DataTree -> SetBranchAddress("mass_mumu", &mass_mumu);
	DataTree -> SetBranchAddress("pt_mumu", &pt_mumu);
	DataTree -> SetBranchAddress("leadMuonP", &leadMuonP);
	DataTree -> SetBranchAddress("leadMuonPt", &leadMuonPt);
	DataTree -> SetBranchAddress("leadMuonEta", &leadMuonEta);
	DataTree -> SetBranchAddress("leadMuonPhi", &leadMuonPhi);
	DataTree -> SetBranchAddress("leadMuonCharge", &leadMuonCharge);
	DataTree -> SetBranchAddress("leadMuonNtrk", &leadMuonNtrk);
	DataTree -> SetBranchAddress("subleadMuonP", &subleadMuonP);
	DataTree -> SetBranchAddress("subleadMuonPt", &subleadMuonPt);
	DataTree -> SetBranchAddress("subleadMuonEta", &subleadMuonEta);
	DataTree -> SetBranchAddress("subleadMuonPhi", &subleadMuonPhi);
	DataTree -> SetBranchAddress("subleadMuonCharge", &subleadMuonCharge);
	DataTree -> SetBranchAddress("subleadMuonNtrk", &subleadMuonNtrk);
	DataTree -> SetBranchAddress("photonE", &photonE);
	DataTree -> SetBranchAddress("photonPt", &photonPt);
	DataTree -> SetBranchAddress("photonEta", &photonEta);
	DataTree -> SetBranchAddress("photonPhi", &photonPhi);
	DataTree -> SetBranchAddress("photonR9", &photonR9);
	DataTree -> SetBranchAddress("photonS4", &photonS4);
	DataTree -> SetBranchAddress("photonEtaWidth", &photonEtaWidth);
	DataTree -> SetBranchAddress("photonSigmaIEtaIEta", &photonSigmaIetaIeta);
	DataTree -> SetBranchAddress("photonEtaEta", &photonSigmaEtaEta);

	nentries = MCTree -> GetEntries();
	RoccoR  rc("rcdata.2016.v3");

	for(int i=0; i<nentries; i++)
	{
		MCTree -> GetEntry(i);
		weight = weight*lumiFactor;
		if(photonPt<30.) continue;

		sumW += weight;
		evMC++;

		if(i%1000==0)
			cout << "Processed " << i << " events out of " << nentries << endl; 

		double mcSFlead = rc.kScaleAndSmearMC(leadMuonCharge, leadMuonPt, leadMuonEta, leadMuonPhi, leadMuonNtrk, gRandom->Rndm(), gRandom->Rndm());
		double mcSFsublead = rc.kScaleAndSmearMC(subleadMuonCharge, subleadMuonPt, subleadMuonEta, subleadMuonPhi, subleadMuonNtrk, gRandom->Rndm(), gRandom->Rndm());

		//double mcSFlead = 1.;
		//double mcSFsublead = 1.;

		mass_mumu = mass_mumu * sqrt (mcSFlead*mcSFsublead);
		mass_mmg = mass_mmg * sqrt (mcSFlead*mcSFsublead);
		pt_mumu = pt_mumu *sqrt (mcSFlead*mcSFsublead);
		leadMuonPt = leadMuonPt * mcSFlead;
		subleadMuonPt = subleadMuonPt * mcSFsublead;
		leadMuonP = leadMuonP * mcSFlead;
		subleadMuonP = subleadMuonP * mcSFsublead;

		MCMass -> Fill(mass_mmg, weight);
		MCMass2 -> Fill(mass_mmg, weight);
		MCMuMuMass -> Fill(mass_mumu, weight);
		MCMuMuPt -> Fill(pt_mumu, weight);
		MCLeadingMuPt -> Fill(leadMuonPt, weight);
		MCLeadingMuEta -> Fill(leadMuonEta, weight);
		MCLeadingMuPhi -> Fill(leadMuonPhi, weight);
		MCLeadingMuCharge -> Fill(leadMuonCharge, weight);
		MCLeadingMuNtrk -> Fill(leadMuonNtrk, weight);
		MCSubleadingMuPt -> Fill(subleadMuonPt, weight);
		MCSubleadingMuEta -> Fill(subleadMuonEta, weight);
		MCSubleadingMuPhi -> Fill(subleadMuonPhi, weight);
		MCSunleadingMuCharge -> Fill(subleadMuonCharge, weight);
		MCSubleadingMuNtrk -> Fill(subleadMuonNtrk, weight);
		MCPhotonE -> Fill(photonE, weight);
		MCPhotonPt -> Fill(photonPt, weight);
		MCPhotonEta -> Fill(photonEta, weight);
		MCPhotonPhi -> Fill(photonPhi, weight);
		MCPhotonR9 -> Fill(photonR9, weight);
		MCPhotonS4 -> Fill(photonS4, weight);
		MCPhotonEtaWidth -> Fill(photonEtaWidth, weight);
		MCPhotonSigmaIetaIeta -> Fill(photonSigmaIetaIeta, weight);
		MCPhotonSigmaEtaEta -> Fill(photonSigmaEtaEta, weight);


		if(fabs(photonEta)<1.4442)
		{
			MCMassEB -> Fill(mass_mmg, weight);
			MCMassEB2 -> Fill(mass_mmg, weight);
			MCPhotonR9EB -> Fill(photonR9, weight);
			MCPhotonS4EB -> Fill(photonS4, weight);
			MCPhotonEtaWidthEB -> Fill(photonEtaWidth, weight);
			MCPhotonSigmaIetaIetaEB -> Fill(photonSigmaIetaIeta, weight);
			MCPhotonSigmaEtaEtaEB -> Fill(photonSigmaEtaEta, weight);

			if(photonR9>0.94)
				MCMassEBHighR9 -> Fill(mass_mmg, weight);
			else
				MCMassEBLowR9 -> Fill(mass_mmg, weight);
		}

		else if(fabs(photonEta)>1.556 && fabs(photonEta)<2.5)
		{
			MCMassEE -> Fill(mass_mmg, weight);
			MCMassEE2 -> Fill(mass_mmg, weight);
			MCPhotonR9EE -> Fill(photonR9, weight);
			MCPhotonS4EE -> Fill(photonS4, weight);
			MCPhotonEtaWidthEE -> Fill(photonEtaWidth, weight);
			MCPhotonSigmaIetaIetaEE -> Fill(photonSigmaIetaIeta, weight);
			MCPhotonSigmaEtaEtaEE -> Fill(photonSigmaEtaEta, weight);

			if(photonR9>0.94)
				MCMassEEHighR9 -> Fill(mass_mmg, weight);
			else
				MCMassEELowR9 -> Fill(mass_mmg, weight);
		}

	}	


	nentries = DataTree -> GetEntries();

	for(int i=0; i<nentries; i++)
	{
		DataTree -> GetEntry(i);
		weight = 1.;

		if(photonPt<30.) continue;

		evData++;
		if(i%1000==0)
			cout << "Processed " << i << " events out of " << nentries << endl;

		double dataSFlead = rc.kScaleDT(leadMuonCharge, leadMuonPt, leadMuonEta, leadMuonPhi); 
		double dataSFsublead = rc.kScaleDT(subleadMuonCharge, subleadMuonPt, subleadMuonEta, subleadMuonPhi);

		//double dataSFlead = 1.; 
		//double dataSFsublead = 1.;

		mass_mumu = mass_mumu * sqrt (dataSFlead*dataSFsublead);
		mass_mmg = mass_mmg * sqrt (dataSFlead*dataSFsublead);
		pt_mumu = pt_mumu *sqrt (dataSFlead*dataSFsublead);
		leadMuonPt = leadMuonPt * dataSFlead;
		subleadMuonPt = subleadMuonPt * dataSFsublead;
		leadMuonP = leadMuonP * dataSFlead;
		subleadMuonP = subleadMuonP * dataSFsublead;

		DataMass -> Fill(mass_mmg, weight);
		DataMass2 -> Fill(mass_mmg, weight);
		DataMuMuMass -> Fill(mass_mumu, weight);
		DataMuMuPt -> Fill(pt_mumu, weight);
		DataLeadingMuPt -> Fill(leadMuonPt, weight);
		DataLeadingMuEta -> Fill(leadMuonEta, weight);
		DataLeadingMuPhi -> Fill(leadMuonPhi, weight);
		DataLeadingMuCharge -> Fill(leadMuonCharge, weight);
		DataLeadingMuNtrk -> Fill(leadMuonNtrk, weight);
		DataSubleadingMuPt -> Fill(subleadMuonPt, weight);
		DataSubleadingMuEta -> Fill(subleadMuonEta, weight);
		DataSubleadingMuPhi -> Fill(subleadMuonPhi, weight);
		DataSunleadingMuCharge -> Fill(subleadMuonCharge, weight);
		DataSubleadingMuNtrk -> Fill(subleadMuonNtrk, weight);
		DataPhotonE -> Fill(photonE, weight);
		DataPhotonPt -> Fill(photonPt, weight);
		DataPhotonEta -> Fill(photonEta, weight);
		DataPhotonPhi -> Fill(photonPhi, weight);
		DataPhotonR9 -> Fill(photonR9, weight);
		DataPhotonS4 -> Fill(photonS4, weight);
		DataPhotonEtaWidth -> Fill(photonEtaWidth, weight);
		DataPhotonSigmaIetaIeta -> Fill(photonSigmaIetaIeta, weight);
		DataPhotonSigmaEtaEta -> Fill(photonSigmaEtaEta, weight);

		MassMuMuVsZ -> Fill(mass_mumu, mass_mmg);
		MassMuMuVsPh -> Fill(mass_mumu, photonPt);
		MassTest -> Fill(mass_mumu, photonPt, mass_mmg);

		if(fabs(photonEta)<1.4442)
		{
			DataMassEB -> Fill(mass_mmg, weight);
			DataMassEB2 -> Fill(mass_mmg, weight);
			DataPhotonR9EB -> Fill(photonR9, weight);
			DataPhotonS4EB -> Fill(photonS4, weight);
			DataPhotonEtaWidthEB -> Fill(photonEtaWidth, weight);
			DataPhotonSigmaIetaIetaEB -> Fill(photonSigmaIetaIeta, weight);
			DataPhotonSigmaEtaEtaEB -> Fill(photonSigmaEtaEta, weight);

			if(photonR9>0.94)
				DataMassEBHighR9 -> Fill(mass_mmg, weight);
			else
				DataMassEBLowR9 -> Fill(mass_mmg, weight);
		}


		else if(fabs(photonEta)>1.556 && fabs(photonEta)<2.5)
		{
			DataMassEE -> Fill(mass_mmg, weight);
			DataMassEE2 -> Fill(mass_mmg, weight);
			DataPhotonR9EE -> Fill(photonR9, weight);
			DataPhotonS4EE -> Fill(photonS4, weight);
			DataPhotonEtaWidthEE -> Fill(photonEtaWidth, weight);
			DataPhotonSigmaIetaIetaEE -> Fill(photonSigmaIetaIeta, weight);
			DataPhotonSigmaEtaEtaEE -> Fill(photonSigmaEtaEta, weight);

			if(photonR9>0.94)
				DataMassEEHighR9 -> Fill(mass_mmg, weight);
			else
				DataMassEELowR9 -> Fill(mass_mmg, weight);
		}

	}	

	for(int binX=1; binX<=MassTest->GetNbinsX(); binX++)
	{
		for(int binY=1; binY<=MassTest->GetNbinsY(); binY++)
		{
			MassTest -> SetBinContent(binX, binY, MassTest-> GetBinContent(binX, binY)/MassMuMuVsPh-> GetBinContent(binX, binY));
		}
	}


	HistoPlotter(MCMass, DataMass, "Mass");
	HistoPlotter(MCMassEB, DataMassEB, "MassEB");
	HistoPlotter(MCMassEE, DataMassEE, "MassEE");
	HistoPlotter(MCMass2, DataMass2, "Mass2");
	HistoPlotter(MCMassEB2, DataMassEB2, "MassEB2");
	HistoPlotter(MCMassEE2, DataMassEE2, "MassEE2");
	HistoPlotter(MCMassEBHighR9, DataMassEBHighR9, "MassEBHighR9");
	HistoPlotter(MCMassEEHighR9, DataMassEEHighR9, "MassEEHighR9");
	HistoPlotter(MCMassEBLowR9, DataMassEBLowR9, "MassEBLowR9");
	HistoPlotter(MCMassEELowR9, DataMassEELowR9, "MassEELowR9");
	HistoPlotter(MCMuMuMass, DataMuMuMass, "MuMuMass");
	HistoPlotter(MCMuMuPt, DataMuMuPt, "MuMuPt");
	HistoPlotter(MCLeadingMuPt, DataLeadingMuPt, "LeadingMuPt");
	HistoPlotter(MCLeadingMuEta, DataLeadingMuEta, "LeadingMuEta");
	HistoPlotter(MCLeadingMuPhi, DataLeadingMuPhi, "LeadingMuPhi");
	HistoPlotter(MCLeadingMuCharge, DataLeadingMuCharge, "LeadingMuCharge");
	HistoPlotter(MCLeadingMuNtrk, DataLeadingMuNtrk, "LeadingMuNtrk");
	HistoPlotter(MCSubleadingMuPt, DataSubleadingMuPt, "SubleadingMuPt");
	HistoPlotter(MCSubleadingMuEta, DataSubleadingMuEta, "SubleadingMuEta");
	HistoPlotter(MCSubleadingMuPhi, DataSubleadingMuPhi, "SubleadingMuPhi");
	HistoPlotter(MCSunleadingMuCharge, DataSunleadingMuCharge, "SunleadingMuCharge");
	HistoPlotter(MCSubleadingMuNtrk, DataSubleadingMuNtrk, "SubleadingMuNtrk");
	HistoPlotter(MCPhotonE, DataPhotonE, "PhotonE");
	HistoPlotter(MCPhotonPt, DataPhotonPt, "PhotonPt");
	HistoPlotter(MCPhotonEta, DataPhotonEta, "PhotonEta");
	HistoPlotter(MCPhotonPhi, DataPhotonPhi, "PhotonPhi");
	HistoPlotter(MCPhotonR9, DataPhotonR9, "PhotonR9");
	HistoPlotter(MCPhotonR9EB, DataPhotonR9EB, "PhotonR9EB");
	HistoPlotter(MCPhotonR9EE, DataPhotonR9EE, "PhotonR9EE");
	HistoPlotter(MCPhotonS4, DataPhotonS4, "PhotonS4");
	HistoPlotter(MCPhotonS4EB, DataPhotonS4EB, "PhotonS4EB");
	HistoPlotter(MCPhotonS4EE, DataPhotonS4EE, "PhotonS4EE");
	HistoPlotter(MCPhotonEtaWidth, DataPhotonEtaWidth, "PhotonEtaWidth");	
	HistoPlotter(MCPhotonEtaWidthEB, DataPhotonEtaWidthEB, "PhotonEtaWidthEB");	
	HistoPlotter(MCPhotonEtaWidthEE, DataPhotonEtaWidthEE, "PhotonEtaWidthEE");	
	HistoPlotter(MCPhotonSigmaIetaIeta, DataPhotonSigmaIetaIeta, "PhotonSigmaIEtaIEta");	
	HistoPlotter(MCPhotonSigmaIetaIetaEB, DataPhotonSigmaIetaIetaEB, "PhotonSigmaIEtaIEtaEB");	
	HistoPlotter(MCPhotonSigmaIetaIetaEE, DataPhotonSigmaIetaIetaEE, "PhotonSigmaIEtaIEtaEE");	
	HistoPlotter(MCPhotonSigmaEtaEta, DataPhotonSigmaEtaEta, "PhotonSigmaEtaEta");	
	HistoPlotter(MCPhotonSigmaEtaEtaEB, DataPhotonSigmaEtaEtaEB, "PhotonSigmaEtaEtaEB");	
	HistoPlotter(MCPhotonSigmaEtaEtaEE, DataPhotonSigmaEtaEtaEE, "PhotonSigmaEtaEtaEE");	
	HistoPlotter2D(MassMuMuVsZ, "MassMuMuVsZ");
	HistoPlotter2D(MassMuMuVsPh, "MassMuMuVsPh");
	HistoPlotter2D(MassTest, "MassTest");


	if(DoNotApplyScaleAndSmearing)
	{
		system("mv *.pdf ~/www/Linearity/Full2016Dataset/MuMuGamma/ControlPlot/NoScale");
		system("mv *.png ~/www/Linearity/Full2016Dataset/MuMuGamma/ControlPlot/NoScale");
	}

	else
	{
		system("mv *.pdf ~/www/Linearity/Full2016Dataset/MuMuGamma/ControlPlot/");
		system("mv *.png ~/www/Linearity/Full2016Dataset/MuMuGamma/ControlPlot/");
	}


	cout << "Data events: " << evData << endl;
	cout << "MC events: " << sumW << " (" << evMC << " unweighted)" << endl; 

}





















