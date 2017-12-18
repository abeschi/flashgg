/*
g++ -Wall -o NewPlotterAfterSelection `root-config --cflags --glibs` -L $ROOTSYS/lib -lRooFitCore -lFoam -lMinuit -lMathMore CMS_lumi.C tdrstyle.C NewPlotterAfterSelection.cpp
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
#include "TAxis.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <cmath>
#include <vector>


using namespace std;

float compute_R(float eta_1, float eta_2, float phi_1, float phi_2)
{
	float delta_eta_squared = (eta_1 - eta_2)*(eta_1 - eta_2);
	float delta_phi = fabs(phi_1 - phi_2);
		if(delta_phi > 3.14159265359)
			delta_phi =2*3.14159265359 - delta_phi;

	return sqrt(delta_eta_squared + delta_phi*delta_phi);
}

float compute_DeltaPhi(float phi_1, float phi_2)
{
	float delta_phi = fabs(phi_1 - phi_2);
		if(delta_phi > 3.14159265359)
			delta_phi =2*3.14159265359 - delta_phi;

	return delta_phi;
}

void MakePlot(TH1F** histos, TString title)
{
	histos[0] -> SetLineWidth(3);			//ttH
	histos[0] -> SetLineColor(kRed + 1);
	histos[0] -> SetFillStyle(0);

	histos[1] -> SetLineWidth(3);			//ggH
	histos[1] -> SetLineColor(kGreen + 2);
	histos[1] -> SetFillStyle(0);

	histos[2] -> SetLineWidth(3);			//VBF
	histos[2] -> SetLineColor(kAzure);
	histos[2] -> SetFillStyle(0);

	histos[3] -> SetLineWidth(3);			//VH
	histos[3] -> SetLineColor(kViolet - 2);
	histos[3] -> SetFillStyle(0);

	histos[4] -> SetLineWidth(3);			//bbH
	histos[4] -> SetLineColor(kOrange);
	histos[4] -> SetFillStyle(0);

	histos[5] -> SetLineWidth(3);			//tHq
	histos[5] -> SetLineColor(kAzure + 8);
	histos[5] -> SetFillStyle(0);

	histos[6] -> SetLineWidth(3);			//tHW
	histos[6] -> SetLineColor(kViolet + 2);
	histos[6] -> SetFillStyle(0);



	histos[7] -> SetMarkerStyle(20);		//Data
	histos[7] -> SetMarkerSize(1);
	histos[7] -> SetMarkerColor(kBlack);
	histos[7] -> SetFillStyle(0);



	histos[8] -> SetLineWidth(1);			//Diphoton
	histos[8] -> SetFillColor(kAzure + 1);
	histos[8] -> SetFillStyle(1001);

	histos[9] -> SetLineWidth(1);			//Gamma + jets
	histos[9] -> SetFillStyle(1001);
	histos[9] -> SetFillColor(kYellow - 4);

	histos[10] -> SetLineWidth(1);			//QCD
	histos[10] -> SetFillColor(kTeal + 9);
	histos[10] -> SetFillStyle(1001);

	histos[11] -> SetLineWidth(1);			//ttGG
	histos[11] -> SetFillColor(kMagenta + 1);
	histos[11] -> SetFillStyle(1001);

	histos[12] -> SetLineWidth(1);			//ttGJets
	histos[12] -> SetFillColor(kMagenta + 1);
	histos[12] -> SetFillStyle(1001);

	histos[13] -> SetLineWidth(1);			//ttJets
	histos[13] -> SetFillColor(kMagenta + 1);
	histos[13] -> SetFillStyle(1001);


	TLegend* leg = new TLegend(0.65, 0.70, 0.9, 0.85);
	leg -> AddEntry(histos[0], "ttH", "l");
	leg -> AddEntry(histos[1], "ggH", "l");
	leg -> AddEntry(histos[2], "VBF", "l");
	leg -> AddEntry(histos[3], "VH", "l");
	leg -> AddEntry(histos[4], "bbH", "l");
	leg -> AddEntry(histos[5], "tHq", "l");
	leg -> AddEntry(histos[6], "tHW", "l");
	leg -> AddEntry(histos[7], "Data sidebands", "p");

	TLegend* leg2 = new TLegend(0.65, 0.70, 0.9, 0.85);
	leg2 -> AddEntry(histos[7], "Data sidebands", "p");
	leg2 -> AddEntry(histos[8], "Diphotons", "f");
	leg2 -> AddEntry(histos[9], "Gamma + jets", "f");
	leg2 -> AddEntry(histos[10], "QCD", "f");
	leg2 -> AddEntry(histos[11], "ttGJets", "f");


	histos[11] -> Add(histos[12]);
	histos[11] -> Add(histos[13]);


	TCanvas* c = new TCanvas();
	c -> cd();

	float m = max(max(histos[0]->GetMaximum()/histos[0]->Integral(), histos[1]->GetMaximum()/histos[1]->Integral()), max(histos[2]->GetMaximum()/histos[2]->Integral(), histos[3]->GetMaximum()/histos[3]->Integral()));
	float m2 = max(max(histos[4]->GetMaximum()/histos[4]->Integral(), histos[5]->GetMaximum()/histos[5]->Integral()), max(histos[6]->GetMaximum()/histos[6]->Integral(), histos[7]->GetMaximum()/histos[7]->Integral()));
	m = max((double)m, (double)m2);

	TH1F* axis = new TH1F(*histos[0]);
	axis -> SetMarkerSize(0);
	axis -> SetLineWidth(0);
	axis -> GetYaxis() -> SetTitleOffset(1.5);
	axis -> GetYaxis() -> SetRangeUser(0, 1.1*m);

	axis -> Draw("histo");
	histos[1] -> DrawNormalized("histo SAME");
	histos[2] -> DrawNormalized("histo SAME");
	histos[3] -> DrawNormalized("histo SAME");
	histos[4] -> DrawNormalized("histo SAME");
	histos[5] -> DrawNormalized("histo SAME");
	histos[6] -> DrawNormalized("histo SAME");
	histos[0] -> DrawNormalized("histo SAME");
	histos[7] -> DrawNormalized("SAME E1");

	leg -> Draw("SAME");
	CMS_lumi(c, 0, 0);

	c -> SaveAs(title + "Signal.png");
	c -> SaveAs(title + "Signal.pdf");


	TCanvas* c2 = new TCanvas();
	c2 -> cd();

	histos[9] -> Add(histos[8]);
	histos[10] -> Add(histos[9]);
	histos[11] -> Add(histos[10]);


	double max2 = histos[7]->GetBinContent(histos[7]->GetMaximumBin());
	if(histos[11]->GetBinContent(histos[11]->GetMaximumBin()) > max2 )
		max2 = histos[11]->GetBinContent(histos[11]->GetMaximumBin());
	histos[11]-> GetYaxis() -> SetRangeUser(0, 1.1*max2);
	histos[11]-> GetYaxis() -> SetTitleSize(0.05);
	histos[11]-> GetYaxis() -> SetTitleFont(42);
	histos[11]-> GetYaxis() -> SetLabelSize(0.045);
	histos[11]-> GetYaxis() -> SetLabelFont(42);
	histos[11]-> GetXaxis() -> SetLabelSize(0);

	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 0.97);
	pad1->SetBottomMargin(0.035);
	pad1->Draw();
	pad1->cd();

	histos[11] -> Draw("histo");
	histos[10] -> Draw("histo SAME");
	histos[9] -> Draw("histo SAME");
	histos[8] -> Draw("histo SAME");
	histos[7] -> Draw("SAME E1");
	leg2 -> Draw("SAME");

	c2->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.10, 1, 0.35);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0);
	pad2->Draw();
	pad2->cd();

	TH1F *h = (TH1F*)histos[7]->Clone("h");
	h->SetLineColor(kBlack);
	h->SetMinimum(0.5);  // Define Y ..
	h->SetMaximum(1.5); // .. range
	h->Sumw2();
	h->SetStats(0);      // No statistics on lower plot
	h->Divide(histos[11]);
	h->SetMarkerStyle(21);
	h -> SetTitle("");
	h-> GetYaxis() -> SetTitle("Data/MC");


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


	CMS_lumi(c2, 0, 10);

	c2 -> SaveAs(title + "Bkg.png");
	c2 -> SaveAs(title + "Bkg.pdf");

	return;
}

int main(int argc, char *argv[])
{
	writeExtraText = true;       // if extra text
	extraText  = "Preliminary";  // default extra text is "Preliminary"
	lumi_sqrtS = "35.9 fb^{-1} (13 TeV)";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

	setTDRStyle();
	gStyle -> SetOptFit(0);
	gStyle -> SetOptStat(0);

	TString SignalFold = "/afs/cern.ch/work/a/abeschi/ttH_MassTop_v5/";
	TString BkgFold = "/afs/cern.ch/work/a/abeschi/ttH_MassTop_Bkg_v3/";
	TString DataFold = "/afs/cern.ch/work/a/abeschi/ttH_MassTop_Data_v2/";
	TString names[14] = {"ttH", "ggH", "VBF", "VH", "bbH", "tHq", "tHW", "Data", "Diphotons", "GJets", "QCD" , "ttGG", "ttGJets", "ttJets"};

	TChain* ttH;
	TChain* ggH;
	TChain* vbf;
	TChain* vh;
	TChain* bbH;
	TChain* tHq;
	TChain* tHW;
	
	TChain* data;

	TChain* diphotons;
	TChain* gjet;
	TChain* qcd;
	TChain* ttGG;
	TChain* ttGJets;
	TChain* ttJets;

	float lumiFactor = 35.9;
	float bDiscriminantThresholdLoose = 0.6;
	float bDiscriminantThresholdMedium = 0.8;
	float bDiscriminantThresholdTight = 0.97;

	bool isLept = (argc!=1 ? 0 : 1);

	if(isLept)
	{	cout << "Processing leptonic tag" << endl;

		ttH = new TChain("TTHLeptonicDumper/trees/tth_13TeV_all");
		ttH -> Add(SignalFold + "ttHJetToGG_M125_13TeV_*.root");		
		ggH = new TChain("TTHLeptonicDumper/trees/ggh_13TeV_all");
		ggH -> Add(SignalFold + "GluGluHToGG_M125_13TeV_*.root");		
		vbf = new TChain("TTHLeptonicDumper/trees/vbf_13TeV_all");
		vbf -> Add(SignalFold + "VBFHToGG_M125_13TeV*.root");		
		vh = new TChain("TTHLeptonicDumper/trees/vh_13TeV_all");
		vh -> Add(SignalFold + "VHToGG_M125_13TeV_*.root");	
		bbH = new TChain("TTHLeptonicDumper/trees/bbh_13TeV_all");
		bbH -> Add(SignalFold + "bbHToGG_M125_*.root");	
		tHq = new TChain("TTHLeptonicDumper/trees/thq_13TeV_all");
		tHq -> Add(SignalFold + "THQ_HToGG_13TeV*.root");	
		tHW = new TChain("TTHLeptonicDumper/trees/thw_13TeV_all");
		tHW -> Add(SignalFold + "THW_HToGG_13TeV*.root");	

		data = new TChain("TTHLeptonicDumper/trees/Data_13TeV_all");
		data -> Add(DataFold + "output_DoubleEG_Run2016*.root");
		
		diphotons = new TChain("TTHLeptonicDumper/trees/diphoton_13TeV_all");
		diphotons -> Add(BkgFold + "output_DiPhotonJetsBox_MGG*.root");
		gjet = new TChain("TTHLeptonicDumper/trees/gjet_13TeV_all");
		gjet -> Add(BkgFold + "output_GJet_Pt*.root");		
		qcd = new TChain("TTHLeptonicDumper/trees/qcd_13TeV_all");
		qcd -> Add(BkgFold + "output_QCD_Pt*.root");
		ttGG = new TChain("TTHLeptonicDumper/trees/ttGG_13TeV_all");
		ttGG -> Add(BkgFold + "output_TTGG_0Jets*.root");
		ttGJets = new TChain("TTHLeptonicDumper/trees/ttGJets_13TeV_all");
		ttGJets -> Add(BkgFold + "output_TTGJets_*.root");
		ttJets = new TChain("TTHLeptonicDumper/trees/ttJets_13TeV_all");
		ttJets -> Add(BkgFold + "output_TTJets_*.root");
	}

	else
	{	cout << "Processing hadronic tag" << endl;

		ttH = new TChain("TTHHadronicDumper/trees/tth_13TeV_all");
		ttH -> Add(SignalFold + "ttHJetToGG_M125_13TeV_*.root");		
		ggH = new TChain("TTHHadronicDumper/trees/ggh_13TeV_all");
		ggH -> Add(SignalFold + "GluGluHToGG_M125_13TeV_*.root");		
		vbf = new TChain("TTHHadronicDumper/trees/vbf_13TeV_all");
		vbf -> Add(SignalFold + "VBFHToGG_M125_13TeV*.root");		
		vh = new TChain("TTHHadronicDumper/trees/vh_13TeV_all");
		vh -> Add(SignalFold + "VHToGG_M125_13TeV_*.root");	
		bbH = new TChain("TTHHadronicDumper/trees/bbh_13TeV_all");
		bbH -> Add(SignalFold + "bbHToGG_M125_*.root");	
		tHq = new TChain("TTHHadronicDumper/trees/thq_13TeV_all");
		tHq -> Add(SignalFold + "THQ_HToGG_13TeV*.root");	
		tHW = new TChain("TTHHadronicDumper/trees/thw_13TeV_all");
		tHW -> Add(SignalFold + "THW_HToGG_13TeV*.root");	

		data = new TChain("TTHHadronicDumper/trees/Data_13TeV_all");
		data -> Add(DataFold + "output_DoubleEG_Run2016*.root");
		
		diphotons = new TChain("TTHHadronicDumper/trees/diphoton_13TeV_all");
		diphotons -> Add(BkgFold + "output_DiPhotonJetsBox_MGG*.root");
		gjet = new TChain("TTHHadronicDumper/trees/gjet_13TeV_all");
		gjet -> Add(BkgFold + "output_GJet_Pt*.root");		
		qcd = new TChain("TTHHadronicDumper/trees/qcd_13TeV_all");
		qcd -> Add(BkgFold + "output_QCD_Pt*.root");
		ttGG = new TChain("TTHHadronicDumper/trees/ttGG_13TeV_all");
		ttGG -> Add(BkgFold + "output_TTGG_0Jets*.root");
		ttGJets = new TChain("TTHHadronicDumper/trees/ttGJets_13TeV_all");
		ttGJets -> Add(BkgFold + "output_TTGJets_*.root");
		ttJets = new TChain("TTHHadronicDumper/trees/ttJets_13TeV_all");
		ttJets -> Add(BkgFold + "output_TTJets_*.root");
	}


	int nvtx = 0;
	float weight = 0;
	
	float dipho_sumpt = 0;
	float dipho_mass = 0;
	float dipho_leadPt = 0;
	float dipho_leadEta = 0;
	float dipho_leadPhi = 0;
	float dipho_leadR9 = 0;
	float dipho_subleadPt = 0;
	float dipho_subleadEta = 0;
	float dipho_subleadPhi = 0;
	float dipho_subleadR9 = 0;
	float dipho_leadIDMVA = 0;
	float dipho_subleadIDMVA = 0;
	float dipho_mva = 0;

	float jet_pt[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
	float jet_eta[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
	float jet_phi[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
	float jet_bdiscriminant[9] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};

	float ele_pt[2] = {0., 0.};
	float ele_eta[2] = {0., 0.};
	float ele_phi[2] = {0., 0.};

	float mu_pt[2] = {0., 0.};
	float mu_eta[2] = {0., 0.};
	float mu_phi[2] = {0., 0.};

	float MetPt = 0.;
	float MetPhi = 0.;
	float ttHMVA = 0.;


	TH1F* nvtx_histo[14];

	TH1F* dipho_sumpt_histo[14];
	TH1F* dipho_mass_histo[14];
	TH1F* dipho_leadPt_histo[14];
	TH1F* dipho_leadPtOverMass_histo[14];
	TH1F* dipho_leadEta_histo[14];
	TH1F* dipho_leadPhi_histo[14];
	TH1F* dipho_leadR9_histo[14];
	TH1F* dipho_subleadPt_histo[14];
	TH1F* dipho_subleadPtOverMass_histo[14];
	TH1F* dipho_subleadEta_histo[14];
	TH1F* dipho_subleadPhi_histo[14];
	TH1F* dipho_subleadR9_histo[14];
	TH1F* dipho_leadIDMVA_histo[14];
	TH1F* dipho_subleadIDMVA_histo[14];
	TH1F* dipho_mva_histo[14];
	TH1F* jet_pt_histo[14];
	TH1F* jet_eta_histo[14];
	TH1F* jet_phi_histo[14];
	TH1F* jet_bdiscriminant_histo[14];
	TH1F* Njet_pt_histo[14][9];
	TH1F* Njet_eta_histo[14][9];
	TH1F* Njet_phi_histo[14][9];
	TH1F* Njet_bdiscriminant_histo[14][9];
	TH1F* DeltaEtaPhotons_histo[14];
	TH1F* DeltaPhiPhotons_histo[14];
	TH1F* DeltaEtaJets01_histo[14];
	TH1F* DeltaPhiJets01_histo[14];
	TH1F* DeltaEtaJets02_histo[14];
	TH1F* DeltaPhiJets02_histo[14];
	TH1F* DeltaEtaJets12_histo[14];
	TH1F* DeltaPhiJets12_histo[14];
	TH1F* njet_histo[14];
	TH1F* nbjet_loose_histo[14];
	TH1F* nbjet_medium_histo[14];
	TH1F* nbjet_tight_histo[14];
	TH1F* nleptons_histo[14];
	TH1F* ele_pt_histo[14];
	TH1F* ele_eta_histo[14];
	TH1F* ele_phi_histo[14];
	TH1F* mu_pt_histo[14];
	TH1F* mu_eta_histo[14];
	TH1F* mu_phi_histo[14];
	TH1F* MetPt_histo[14];
	TH1F* MetPhi_histo[14];
	TH1F* ttHMVA_histo[14];

	for(int i=0; i<14; i++)
	{	
		nvtx_histo[i] = new TH1F(("nvtx_histo"+ std::to_string(i)).c_str(), "; n_{vtx}; Counts", 60, -0.5, 59.5 );
		
		dipho_sumpt_histo[i] = new TH1F(("dipho_sumpt_histo"+ std::to_string(i)).c_str(), "; P_{T}^{#gamma#gamma} (GeV); Counts", 160, 80, 300 );
		dipho_mass_histo[i] = new TH1F(("dipho_mass_histo"+ std::to_string(i)).c_str(), "; m_{#gamma#gamma} (GeV); Counts", 80, 100, 180 );
		dipho_leadPt_histo[i] = new TH1F(("dipho_leadPt_histo"+ std::to_string(i)).c_str(), "; P_{T}^{leading #gamma} (GeV); Counts", 90, 20, 200 );
		dipho_leadPtOverMass_histo[i] = new TH1F(("dipho_leadPtOverMass_histo"+ std::to_string(i)).c_str(), "; P_{T}^{leading #gamma}/mass_{#gamma #gamma} (GeV); Counts", 40, 0, 1.5 );
		dipho_leadEta_histo[i] = new TH1F(("dipho_leadEta_histo"+ std::to_string(i)).c_str(), "; #eta; Counts", 50, -2.5, 2.5 );
		dipho_leadPhi_histo[i] = new TH1F(("dipho_leadPhi_histo"+ std::to_string(i)).c_str(), "; #varphi; Counts", 50, -3.15, 3.15 );
		dipho_leadR9_histo[i] = new TH1F(("dipho_leadR9_histo"+ std::to_string(i)).c_str(), "; R_{9}; Counts", 50, 0.7, 1 );
		dipho_subleadPt_histo[i] = new TH1F(("dipho_subleadPt_histo"+ std::to_string(i)).c_str(), "; P_{T}^{subleading #gamma} (GeV); Counts", 65, 20, 150 );
		dipho_subleadPtOverMass_histo[i] = new TH1F(("dipho_subleadPtOverMass_histo"+ std::to_string(i)).c_str(), "; P_{T}^{subleading #gamma}/mass_{#gamma #gamma} (GeV); Counts", 50, 0, 2 );
		dipho_subleadEta_histo[i] = new TH1F(("dipho_subleadEta_histo"+ std::to_string(i)).c_str(), "; #eta; Counts", 50, -2.5, 2.5 );
		dipho_subleadPhi_histo[i] = new TH1F(("dipho_subleadPhi_histo"+ std::to_string(i)).c_str(), "; #varphi; Counts", 50, -3.15, 3.15 );
		dipho_subleadR9_histo[i] = new TH1F(("dipho_subleadR9_histo"+ std::to_string(i)).c_str(), "; R_{9}; Counts", 50, 0.7, 1 );
		dipho_leadIDMVA_histo[i] = new TH1F(("dipho_leadIDMVA_histo"+ std::to_string(i)).c_str(), "; Photon IDMVA; Counts", 50, -1, 1 );
		dipho_subleadIDMVA_histo[i] = new TH1F(("dipho_subleadIDMVA_histo"+ std::to_string(i)).c_str(), "; Photon IDMVA; Counts", 50, -1, 1 );
		dipho_mva_histo[i] = new TH1F(("dipho_mva_histo"+ std::to_string(i)).c_str(), "; Diphoton MVA; Counts", 50, -1, 1 );

		jet_pt_histo[i] = new TH1F(("jet_pt_histo"+ std::to_string(i)).c_str(), "; Jet p_{T} (GeV); Counts", 100, 0, 200 );
		jet_eta_histo[i] = new TH1F(("jet_eta_histo"+ std::to_string(i)).c_str(), "; Jet #eta; Counts", 50, -2.5, 2.5 );
		jet_phi_histo[i] = new TH1F(("jet_phi_histo"+ std::to_string(i)).c_str(), "; Jet #varphi; Counts", 50, -3.15, 3.15 );
		jet_bdiscriminant_histo[i] = new TH1F(("jet_bdiscriminant_histo"+ std::to_string(i)).c_str(), "; Jet bdiscriminant; Counts", 50, 0, 1 );

		DeltaEtaPhotons_histo[i] = new TH1F(("DeltaEtaPhotons_histo"+ std::to_string(i)).c_str(), "; #Delta#eta^{#gamma #gamma}; Counts", 50, -3, 3 );
		DeltaPhiPhotons_histo[i] = new TH1F(("DeltaPhiPhotons_histo"+ std::to_string(i)).c_str(), "; #Delta#varphi^{#gamma #gamma}; Counts", 50, 0, 3.15 );
		DeltaEtaJets01_histo[i] = new TH1F(("DeltaEtaJets01_histo"+ std::to_string(i)).c_str(), "; #Delta#eta^{Jet0 Jet1}; Counts", 50, -3, 3 );
		DeltaPhiJets01_histo[i] = new TH1F(("DeltaPhiJets01_histo"+ std::to_string(i)).c_str(), "; #Delta#varphi^{Jet0 Jet1}; Counts", 50, 0, 3.15 );
		DeltaEtaJets02_histo[i] = new TH1F(("DeltaEtaJets02_histo"+ std::to_string(i)).c_str(), "; #Delta#eta^{Jet0 Jet2}; Counts", 50, -3, 3 );
		DeltaPhiJets02_histo[i] = new TH1F(("DeltaPhiJets02_histo"+ std::to_string(i)).c_str(), "; #Delta#varphi^{Jet0 Jet2}; Counts", 50, 0, 3.15 );
		DeltaEtaJets12_histo[i] = new TH1F(("DeltaEtaJets12_histo"+ std::to_string(i)).c_str(), "; #Delta#eta^{Jet1 Jet2}; Counts", 50, -3, 3 );
		DeltaPhiJets12_histo[i] = new TH1F(("DeltaPhiJets12_histo"+ std::to_string(i)).c_str(), "; #Delta#varphi^{Jet1 Jet2}; Counts", 50, 0, 3.15 );
		
		for(int j=0; j<9; j++)	
		{
			Njet_pt_histo[i][j] = new TH1F(("leadingjet_pt_histo"+ std::to_string(i) + std::to_string(j)).c_str(), "; Jet p_{T} (GeV); Counts", 100, 0, 200 );
			Njet_eta_histo[i][j] = new TH1F(("leadingjet_eta_histo"+ std::to_string(i) + std::to_string(j)).c_str(), "; Jet #eta; Counts", 50, -2.5, 2.5 );
			Njet_phi_histo[i][j] = new TH1F(("leadingjet_phi_histo"+ std::to_string(i) + std::to_string(j)).c_str(), "; Jet #varphi; Counts", 50, -3.15, 3.15 );
			Njet_bdiscriminant_histo[i][j] = new TH1F(("leadingjet_bdiscriminant_histo"+ std::to_string(i) + std::to_string(j)).c_str(), "; Jet bdiscriminant; Counts", 50, 0, 1 );
		}

		njet_histo[i] = new TH1F(("njet_histo"+ std::to_string(i)).c_str(), "; Number of jets; Counts", 16, -0.5, 15.5);
		nbjet_loose_histo[i] = new TH1F(("nbjet_loose_histo"+ std::to_string(i)).c_str(), "; Number of loose b-jets; Counts", 6, -0.5, 5.5);
		nbjet_medium_histo[i] = new TH1F(("nbjet_medium_histo"+ std::to_string(i)).c_str(), "; Number of medium b-jets; Counts", 6, -0.5, 5.5);
		nbjet_tight_histo[i] = new TH1F(("nbjet_tight_histo"+ std::to_string(i)).c_str(), "; Number of tight b-jets; Counts", 6, -0.5, 5.5);
		nleptons_histo[i] = new TH1F(("nleptons_histo"+ std::to_string(i)).c_str(), "; Number of leptons; Counts", 6, -0.5, 5.5);

		ele_pt_histo[i] = new TH1F(("ele_pt_histo"+ std::to_string(i)).c_str(), "; Ele p_{T} (GeV); Counts", 200, 0, 200 );
		ele_eta_histo[i] = new TH1F(("ele_eta_histo"+ std::to_string(i)).c_str(), "; Ele #eta; Counts", 100, -2.5, 2.5 );
		ele_phi_histo[i] = new TH1F(("ele_phi_histo"+ std::to_string(i)).c_str(), "; Ele #varphi; Counts", 100, -3.15, 3.15 );
		mu_pt_histo[i] = new TH1F(("mu_pt_histo"+ std::to_string(i)).c_str(), "; Mu p_{T} (GeV); Counts", 200, 0, 200 );
		mu_eta_histo[i] = new TH1F(("mu_eta_histo"+ std::to_string(i)).c_str(), "; Mu #eta; Counts", 100, -2.5, 2.5 );
		mu_phi_histo[i] = new TH1F(("mu_phi_histo"+  std::to_string(i)).c_str(), "; Mu #varphi; Counts", 100, -3.15, 3.15 );

		MetPt_histo[i] = new TH1F(("MetPt_histo"+  std::to_string(i)).c_str(), "; P_{T}^{MET} (GeV); Counts", 100, 0, 200 );
		MetPhi_histo[i] = new TH1F(("MetPhi_histo"+  std::to_string(i)).c_str(), "; #varphi^{MET}; Counts", 50, -3.15, 3.15 );
		ttHMVA_histo[i] = new TH1F(("ttHMVA_histo"+  std::to_string(i)).c_str(), "; ttHMVA; Counts", 50, -1, 1 );
	}

	for(int n=0; n<14; n++)
	{
		int nentries;
		TChain* serviceTree;
		switch (n)
		{
			case(0):
				nentries = ttH -> GetEntries();
				serviceTree = (TChain*)ttH -> Clone();
				break;


			case(1):
				nentries = ggH -> GetEntries();
				serviceTree = (TChain*)ggH -> Clone();
				break;


			case(2):
				nentries = vbf -> GetEntries();
				serviceTree = (TChain*)vbf -> Clone();
				break;


			case(3):
				nentries = vh -> GetEntries();
				serviceTree = (TChain*)vh -> Clone();
				break;

			case(4):
				nentries = bbH -> GetEntries();
				serviceTree = (TChain*)bbH -> Clone();
				break;

			case(5):
				nentries = tHq -> GetEntries();
				serviceTree = (TChain*)tHq -> Clone();
				break;

			case(6):
				nentries = tHW -> GetEntries();
				serviceTree = (TChain*)tHW -> Clone();
				break;

			case(7):
				nentries = data -> GetEntries();
				serviceTree = (TChain*)data -> Clone();
				break;


			case(8):
				nentries = diphotons -> GetEntries();
				serviceTree = (TChain*)diphotons -> Clone();
				break;

			case(9):
				nentries = gjet -> GetEntries();
				serviceTree = (TChain*)gjet -> Clone();
				break;


			case(10):
				nentries = qcd -> GetEntries();
				serviceTree = (TChain*)qcd -> Clone();
				break;

			case(11):
				nentries = ttGG -> GetEntries();
				serviceTree = (TChain*)ttGG -> Clone();
				break;

			case(12):
				nentries = ttGJets -> GetEntries();
				serviceTree = (TChain*)ttGJets -> Clone();
				break;

			case(13):
				nentries = ttJets -> GetEntries();
				serviceTree = (TChain*)ttJets -> Clone();
				break;

			default:
				nentries = 0;
		}

		serviceTree -> SetBranchAddress("nvtx", &nvtx);
		serviceTree -> SetBranchAddress("weight", &weight);
		serviceTree -> SetBranchAddress("dipho_sumpt", &dipho_sumpt);
		serviceTree -> SetBranchAddress("dipho_mass", &dipho_mass);
		serviceTree -> SetBranchAddress("dipho_leadPt", &dipho_leadPt);
		serviceTree -> SetBranchAddress("dipho_leadEta", &dipho_leadEta);
		serviceTree -> SetBranchAddress("dipho_leadPhi", &dipho_leadPhi);
		serviceTree -> SetBranchAddress("dipho_leadR9", &dipho_leadR9);
		serviceTree -> SetBranchAddress("dipho_subleadPt", &dipho_subleadPt);
		serviceTree -> SetBranchAddress("dipho_subleadEta", &dipho_subleadEta);
		serviceTree -> SetBranchAddress("dipho_subleadPhi", &dipho_subleadPhi);
		serviceTree -> SetBranchAddress("dipho_subleadR9", &dipho_subleadR9);
		serviceTree -> SetBranchAddress("dipho_leadIDMVA", &dipho_leadIDMVA);
		serviceTree -> SetBranchAddress("dipho_subleadIDMVA", &dipho_subleadIDMVA);
		serviceTree -> SetBranchAddress("dipho_mva", &dipho_mva);


		if(isLept)
		{
			for(int i=1; i<6; i++)
			{	serviceTree -> SetBranchAddress(("jet_pt"+ std::to_string(i)).c_str(), &jet_pt[i-1]);
				serviceTree -> SetBranchAddress(("jet_eta"+ std::to_string(i)).c_str(), &jet_eta[i-1]);
				serviceTree -> SetBranchAddress(("jet_phi"+ std::to_string(i)).c_str(), &jet_phi[i-1]);
				serviceTree -> SetBranchAddress(("jet_bdiscriminant"+ std::to_string(i)).c_str(), &jet_bdiscriminant[i-1]);

			}

			for(int i=1; i<3; i++)
			{	serviceTree -> SetBranchAddress(("ele_pt"+ std::to_string(i)).c_str(), &ele_pt[i-1]);
				serviceTree -> SetBranchAddress(("ele_eta"+ std::to_string(i)).c_str(), &ele_eta[i-1]);
				serviceTree -> SetBranchAddress(("ele_phi"+ std::to_string(i)).c_str(), &ele_phi[i-1]);
				serviceTree -> SetBranchAddress(("mu_pt"+ std::to_string(i)).c_str(), &mu_pt[i-1]);
				serviceTree -> SetBranchAddress(("mu_eta"+ std::to_string(i)).c_str(), &mu_eta[i-1]);
				serviceTree -> SetBranchAddress(("mu_phi"+ std::to_string(i)).c_str(), &mu_phi[i-1]);
			}
			serviceTree -> SetBranchAddress("MetPt", &MetPt);
			serviceTree -> SetBranchAddress("MetPhi", &MetPhi);
		}
		else
		{
			for(int i=1; i<10; i++)
			{	serviceTree -> SetBranchAddress(("jet_pt"+ std::to_string(i)).c_str(), &jet_pt[i-1]);
				serviceTree -> SetBranchAddress(("jet_eta"+ std::to_string(i)).c_str(), &jet_eta[i-1]);
				serviceTree -> SetBranchAddress(("jet_phi"+ std::to_string(i)).c_str(), &jet_phi[i-1]);

				serviceTree -> SetBranchAddress(("jet_bdiscriminant"+ std::to_string(i)).c_str(), &jet_bdiscriminant[i-1]);
			}

			serviceTree -> SetBranchAddress("ttHMVA", &ttHMVA);
			serviceTree -> SetBranchAddress("MetPt", &MetPt);
			serviceTree -> SetBranchAddress("MetPhi", &MetPhi);
		}


		for(int i=0; i<nentries; i++)
//		for(int i=0; i<100; i++)
		{	
			serviceTree -> GetEntry(i);
			if(i%10000==0) cout << "Processing tag " << names[n] << ", event " << i << " out of " << nentries << "\r" << flush;

			int nleptons = 0;
			int njet = 0;
			int nbjet_loose = 0; 
			int nbjet_medium = 0; 
			int nbjet_tight = 0; 

			if(dipho_mass<100 || dipho_mass>180) continue;
			if(!isLept && ttHMVA<0.75) continue;
			if(n==7 && (dipho_mass>115 && dipho_mass<135)) continue;
			if(n==7)
				weight = 1./lumiFactor;

			dipho_mass_histo[n] -> Fill(dipho_mass, weight*lumiFactor);
			if(n>7 && (dipho_mass>115 && dipho_mass<135)) continue;

			nvtx_histo[n] -> Fill(nvtx, weight*lumiFactor);
			dipho_sumpt_histo[n] -> Fill(dipho_sumpt, weight*lumiFactor);
			dipho_leadPt_histo[n] -> Fill(dipho_leadPt, weight*lumiFactor);
			dipho_leadPtOverMass_histo[n] -> Fill(dipho_leadPt/dipho_mass, weight*lumiFactor);
			dipho_leadEta_histo[n] -> Fill(dipho_leadEta, weight*lumiFactor);
			dipho_leadPhi_histo[n] -> Fill(dipho_leadPhi, weight*lumiFactor);
			dipho_leadR9_histo[n] -> Fill(dipho_leadR9, weight*lumiFactor);
			dipho_subleadPt_histo[n] -> Fill(dipho_subleadPt, weight*lumiFactor);
			dipho_subleadPtOverMass_histo[n] -> Fill(dipho_subleadPt/dipho_mass, weight*lumiFactor);
			dipho_subleadEta_histo[n] -> Fill(dipho_subleadEta, weight*lumiFactor);
			dipho_subleadPhi_histo[n] -> Fill(dipho_subleadPhi, weight*lumiFactor);
			dipho_subleadR9_histo[n] -> Fill(dipho_subleadR9, weight*lumiFactor);
			dipho_leadIDMVA_histo[n] -> Fill(dipho_leadIDMVA, weight*lumiFactor);
			dipho_subleadIDMVA_histo[n] -> Fill(dipho_subleadIDMVA, weight*lumiFactor);
			dipho_mva_histo[n] -> Fill(dipho_mva, weight*lumiFactor);
			MetPt_histo[n] -> Fill(MetPt, weight*lumiFactor);
			MetPhi_histo[n] -> Fill(MetPhi, weight*lumiFactor);
			DeltaEtaPhotons_histo[n] -> Fill(dipho_leadEta - dipho_subleadEta, weight*lumiFactor);
			DeltaPhiPhotons_histo[n] -> Fill(compute_DeltaPhi(dipho_leadPhi, dipho_subleadPhi), weight*lumiFactor);
                        DeltaEtaJets01_histo[n] -> Fill(jet_eta[0] - jet_eta[1], weight*lumiFactor);
                        DeltaPhiJets01_histo[n] -> Fill(compute_DeltaPhi(jet_phi[0], jet_phi[1]), weight*lumiFactor);
			if(jet_pt[2]>0.)
			{
                                DeltaEtaJets02_histo[n] -> Fill(jet_eta[0] - jet_eta[2], weight*lumiFactor);
                                DeltaPhiJets02_histo[n] -> Fill(compute_DeltaPhi(jet_phi[0], jet_phi[2]), weight*lumiFactor);
                                DeltaEtaJets12_histo[n] -> Fill(jet_eta[1] - jet_eta[2], weight*lumiFactor);
                                DeltaPhiJets12_histo[n] -> Fill(compute_DeltaPhi(jet_phi[1], jet_phi[2]), weight*lumiFactor);
                        }

			if(isLept)
			{
				for(int jIndex=0; jIndex<5; jIndex++)
				{
					if(jet_pt[jIndex]>0.)
					{
						jet_pt_histo[n] -> Fill(jet_pt[jIndex], weight*lumiFactor);
						jet_eta_histo[n] -> Fill(jet_eta[jIndex], weight*lumiFactor);
						jet_phi_histo[n] -> Fill(jet_phi[jIndex], weight*lumiFactor);
						jet_bdiscriminant_histo[n] -> Fill(jet_bdiscriminant[jIndex], weight*lumiFactor);
						Njet_pt_histo[n][jIndex] -> Fill(jet_pt[jIndex], weight*lumiFactor);
						Njet_eta_histo[n][jIndex] -> Fill(jet_eta[jIndex], weight*lumiFactor);
						Njet_phi_histo[n][jIndex] -> Fill(jet_phi[jIndex], weight*lumiFactor);
						Njet_bdiscriminant_histo[n][jIndex] -> Fill(jet_bdiscriminant[jIndex], weight*lumiFactor);


						njet++;
						if(jet_bdiscriminant[jIndex]>bDiscriminantThresholdLoose)
							nbjet_loose++;
						if(jet_bdiscriminant[jIndex]>bDiscriminantThresholdMedium)
							nbjet_medium++;
						if(jet_bdiscriminant[jIndex]>bDiscriminantThresholdTight)
							nbjet_tight++;
					}
				}

				if(ele_pt[0]>0.)
				{
					ele_pt_histo[n] -> Fill(ele_pt[0], weight*lumiFactor);
					ele_eta_histo[n] -> Fill(ele_eta[0], weight*lumiFactor);
					ele_phi_histo[n] -> Fill(ele_phi[0], weight*lumiFactor);

					nleptons++;
				}

				if(ele_pt[1]>0.)
				{
					ele_pt_histo[n] -> Fill(ele_pt[1], weight*lumiFactor);
					ele_eta_histo[n] -> Fill(ele_eta[1], weight*lumiFactor);
					ele_phi_histo[n] -> Fill(ele_phi[1], weight*lumiFactor);

					nleptons++;
				}

				if(mu_pt[0]>0.)
				{
					mu_pt_histo[n] -> Fill(mu_pt[0], weight*lumiFactor);
					mu_eta_histo[n] -> Fill(mu_eta[0], weight*lumiFactor);
					mu_phi_histo[n] -> Fill(mu_phi[0], weight*lumiFactor);

					nleptons++;
				}

				if(mu_pt[1]>0.)
				{
					mu_pt_histo[n] -> Fill(mu_pt[1], weight*lumiFactor);
					mu_eta_histo[n] -> Fill(mu_eta[1], weight*lumiFactor);
					mu_phi_histo[n] -> Fill(mu_phi[1], weight*lumiFactor);

					nleptons++;
				}
			}

			else
			{
				for(int jIndex=0; jIndex<9; jIndex++)
				{
					if(jet_pt[jIndex]>0.)
					{
						jet_pt_histo[n] -> Fill(jet_pt[jIndex], weight*lumiFactor);
						jet_eta_histo[n] -> Fill(jet_eta[jIndex], weight*lumiFactor);
						jet_phi_histo[n] -> Fill(jet_phi[jIndex], weight*lumiFactor);
						jet_bdiscriminant_histo[n] -> Fill(jet_bdiscriminant[jIndex], weight*lumiFactor);
						Njet_pt_histo[n][jIndex] -> Fill(jet_pt[jIndex], weight*lumiFactor);
						Njet_eta_histo[n][jIndex] -> Fill(jet_eta[jIndex], weight*lumiFactor);
						Njet_phi_histo[n][jIndex] -> Fill(jet_phi[jIndex], weight*lumiFactor);
						Njet_bdiscriminant_histo[n][jIndex] -> Fill(jet_bdiscriminant[jIndex], weight*lumiFactor);

						njet++;
						if(jet_bdiscriminant[jIndex]>bDiscriminantThresholdLoose)
							nbjet_loose++;
						if(jet_bdiscriminant[jIndex]>bDiscriminantThresholdMedium)
							nbjet_medium++;
						if(jet_bdiscriminant[jIndex]>bDiscriminantThresholdTight)
							nbjet_tight++;
					}
				}

				ttHMVA_histo[n] -> Fill(ttHMVA, weight*lumiFactor);
			}


			njet_histo[n] -> Fill(njet, weight*lumiFactor);
			nbjet_loose_histo[n] -> Fill(nbjet_loose, weight*lumiFactor);
			nbjet_medium_histo[n] -> Fill(nbjet_medium, weight*lumiFactor);
			nbjet_tight_histo[n] -> Fill(nbjet_tight, weight*lumiFactor);
			nleptons_histo[n] -> Fill(nleptons, weight*lumiFactor);
		}

		delete serviceTree;
		cout << "Processed tag " << names[n] << ", " << nentries << " events out of " << nentries << endl;

	}

	MakePlot(nvtx_histo, "Vertex");
	MakePlot(dipho_sumpt_histo, "DiphotonSumpt");
	MakePlot(dipho_mass_histo, "DiphotonMass");
	MakePlot(dipho_leadPt_histo, "LeadingPhotonPt");
	MakePlot(dipho_leadPtOverMass_histo, "LeadingPhotonPtOverDiphotonMass");
	MakePlot(dipho_leadEta_histo, "LeadingPhotonEta");
	MakePlot(dipho_leadPhi_histo, "LeadingPhotonPhi");
	MakePlot(dipho_leadR9_histo, "LeadingPhotonR9");
	MakePlot(dipho_subleadPt_histo, "SubleadingPhotonPt");
	MakePlot(dipho_subleadPtOverMass_histo, "SubleadingPhotonPtOverDiphotnMass");
	MakePlot(dipho_subleadEta_histo, "SubleadingPhotonEta");
	MakePlot(dipho_subleadPhi_histo, "SubleadingPhotonPhi");
	MakePlot(dipho_subleadR9_histo, "SubleadingPhotonR9");
	MakePlot(dipho_leadIDMVA_histo, "LeadingPhotonIDMVA");
	MakePlot(dipho_subleadIDMVA_histo, "SubleadingPhotonIDMVA");
	MakePlot(dipho_mva_histo, "DiphotonMVA");
	MakePlot(jet_pt_histo, "JetPt");
	MakePlot(jet_eta_histo, "JetEta");
	MakePlot(jet_phi_histo, "JetPhi");
	MakePlot(jet_bdiscriminant_histo, "JetBdiscrimiant");
	MakePlot(MetPt_histo, "MetPt");
	MakePlot(MetPhi_histo, "MetPhi");
	MakePlot(njet_histo, "Njet");
	MakePlot(nbjet_loose_histo, "NbjetLoose");
	MakePlot(nbjet_medium_histo, "NbjetMedium");
	MakePlot(nbjet_tight_histo, "NbjetTight");
	MakePlot(DeltaEtaPhotons_histo, "DeltaEtaPhotons");
	MakePlot(DeltaPhiPhotons_histo, "DeltaPhiPhotons");
	MakePlot(DeltaEtaJets01_histo, "DeltaEtaJets01");
	MakePlot(DeltaPhiJets01_histo, "DeltaPhiJets01");
	MakePlot(DeltaEtaJets02_histo, "DeltaEtaJets02");
	MakePlot(DeltaPhiJets02_histo, "DeltaPhiJets02");
	MakePlot(DeltaEtaJets12_histo, "DeltaEtaJets12");
	MakePlot(DeltaPhiJets12_histo, "DeltaPhiJets12");

	if(isLept)
	{	MakePlot(ele_pt_histo, "ElePt");
		MakePlot(ele_eta_histo, "EleEta");
		MakePlot(ele_phi_histo, "ElePhi");
		MakePlot(mu_pt_histo, "MuPt");
		MakePlot(mu_eta_histo, "MuEta");
		MakePlot(mu_phi_histo, "MuPhi");
		MakePlot(nleptons_histo, "Nleptons");
	}

	else
	{
		MakePlot(ttHMVA_histo, "ttHMVA");
	}



	if(isLept)
	{	system("mv *.png ~/www/ttH/Full2016Dataset/Leptonic/AfterSelection");
		system("mv *.pdf ~/www/ttH/Full2016Dataset/Leptonic/AfterSelection");
	}
	else
	{	system("mv *.png ~/www/ttH/Full2016Dataset/Hadronic/AfterSelection");
		system("mv *.pdf ~/www/ttH/Full2016Dataset/Hadronic/AfterSelection");
	}



	if(isLept)
	{	
		for(int i=0; i<5; i++)
		{
			TH1F* tmp_pt[14];
			TH1F* tmp_eta[14];
			TH1F* tmp_phi[14];
			TH1F* tmp_bdiscriminant[14];

			for(int j=0; j<14; j++)
			{	tmp_pt[j] = Njet_pt_histo[j][i];
				tmp_eta[j] = Njet_eta_histo[j][i];
				tmp_phi[j] = Njet_phi_histo[j][i];
				tmp_bdiscriminant[j] = Njet_bdiscriminant_histo[j][i];
			}

			MakePlot(tmp_pt, ("JetPt"+ to_string(i)).c_str());
			MakePlot(tmp_eta, ("JetEta"+ to_string(i)).c_str());
			MakePlot(tmp_phi, ("JetPhi"+ to_string(i)).c_str());
			MakePlot(tmp_bdiscriminant, ("JetBdiscriminant"+ to_string(i)).c_str());
		}
	}

	else
	{
		for(int i=0; i<9; i++)
		{
			TH1F* tmp_pt[14];
			TH1F* tmp_eta[14];
			TH1F* tmp_phi[14];
			TH1F* tmp_bdiscriminant[14];

			for(int j=0; j<14; j++)
			{	tmp_pt[j] = Njet_pt_histo[j][i];
				tmp_eta[j] = Njet_eta_histo[j][i];
				tmp_phi[j] = Njet_phi_histo[j][i];
				tmp_bdiscriminant[j] = Njet_bdiscriminant_histo[j][i];
			}

			MakePlot(tmp_pt, ("JetPt"+ to_string(i)).c_str());
			MakePlot(tmp_eta, ("JetEta"+ to_string(i)).c_str());
			MakePlot(tmp_phi, ("JetPhi"+ to_string(i)).c_str());
			MakePlot(tmp_bdiscriminant, ("JetBdiscriminant"+ to_string(i)).c_str());
		}

	}


	if(isLept)
	{	system("mv *.png ~/www/ttH/Full2016Dataset/Leptonic/AfterSelection/Jets");
		system("mv *.pdf ~/www/ttH/Full2016Dataset/Leptonic/AfterSelection/Jets");
	}
	else
	{	system("mv *.png ~/www/ttH/Full2016Dataset/Hadronic/AfterSelection/Jets");
		system("mv *.pdf ~/www/ttH/Full2016Dataset/Hadronic/AfterSelection/Jets");
	}


	cout << endl << "Lufthansa, partner of Star Alliance, thanks you for choosign our company, we hope to see you again on board of our airctafts" << endl << endl;



























}
