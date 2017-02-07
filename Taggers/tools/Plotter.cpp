/*
g++ -Wall -o Plotter `root-config --cflags --glibs` -L $ROOTSYS/lib -lRooFitCore -lFoam -lMinuit -lMathMore CMS_lumi.C tdrstyle.C Plotter.cpp
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

float NormalizationFactor(TH1F* h)
{

	float n = 1;

	TF1* f = new TF1("f", "expo", 100, 180);

	int i=0;
	int counter = 0;

	while(counter<10 && i!=0)
	{
		TFitResultPtr p = h -> Fit("f", "RQN");
		i = p;
		counter++;
	}

	
	if(i==0)
		n = 1. - f->Integral(115, 135)/f->Integral(100, 180);

	return n;
}




void MakePlot(TH1F** histos, TString title, int Pos=0)
{
	histos[0] -> SetLineWidth(3);
	histos[0] -> SetLineColor(kRed);
	histos[0] -> SetFillStyle(0);

	histos[1] -> SetLineWidth(3);
	histos[1] -> SetLineColor(kGreen + 2);
	histos[1] -> SetFillStyle(0);

	histos[2] -> SetLineWidth(3);
	histos[2] -> SetLineColor(kAzure);
	histos[2] -> SetFillStyle(0);

	histos[3] -> SetLineWidth(3);
	histos[3] -> SetLineColor(kViolet-2);
	histos[3] -> SetFillStyle(0);

	histos[4] -> SetMarkerStyle(20);
	histos[4] -> SetMarkerSize(1);
	histos[4] -> SetMarkerColor(kBlack);
	histos[4] -> SetFillStyle(0);

	histos[5] -> SetLineWidth(1);
	histos[5] -> SetFillColor(kAzure + 1);
	histos[5] -> SetFillStyle(1001);

	histos[6] -> SetLineWidth(1);
	histos[6] -> SetFillColor(kYellow - 4);
	histos[6] -> SetFillStyle(1001);

	histos[7] -> SetLineWidth(1);
	histos[7] -> SetFillColor(kTeal + 9);
	histos[7] -> SetFillStyle(1001);



	TLegend* leg = new TLegend(0.73, 0.73, 0.92, 0.9);
	leg -> AddEntry(histos[0], "ttH", "l");
	leg -> AddEntry(histos[1], "ggH", "l");
	leg -> AddEntry(histos[2], "vbf", "l");
	leg -> AddEntry(histos[3], "vh", "l");
	leg -> AddEntry(histos[4], "Data sidebands", "p");

	TLegend* leg2 = new TLegend(0.73, 0.70, 0.92, 0.9);
	leg2 -> AddEntry(histos[4], "Data sidebands", "p");
	leg2 -> AddEntry(histos[5], "Diphotons", "f");
	leg2 -> AddEntry(histos[6], "Gamma + jets", "f");
	leg2 -> AddEntry(histos[7], "QCD", "f");


	if(Pos==1)
	{
		leg -> SetX1(0.23);
		leg -> SetX2(0.4);
		leg2 -> SetX1(0.23);
		leg2 -> SetX2(0.4);
		leg -> SetY1(0.63);
		leg -> SetY2(0.8);
		leg2 -> SetY1(0.58);
		leg2 -> SetY2(0.77);

	}

	if(Pos==2)
	{
		leg -> SetX1(0.23);
		leg -> SetX2(0.42);
		leg2 -> SetX1(0.23);
		leg2 -> SetX2(0.42);
		leg -> SetY1(0.23);
		leg -> SetY2(0.4);
		leg2 -> SetY1(0.23);
		leg2 -> SetY2(0.4);

	}

	TCanvas* c = new TCanvas();
	c -> cd();

	float m = max(max(histos[0]->GetMaximum()/histos[0]->Integral(), histos[1]->GetMaximum()/histos[1]->Integral()), max(histos[2]->GetMaximum()/histos[2]->Integral(), histos[3]->GetMaximum()/histos[3]->Integral()));
	m = max((double)m, histos[4]->GetMaximum()/histos[4]->Integral());


	TH1F* axis = new TH1F(*histos[0]);
	axis -> SetMarkerSize(0);
	axis -> SetLineWidth(0);
	axis -> GetYaxis() -> SetTitleOffset(1.5);
	axis -> GetYaxis() -> SetRangeUser(0, 1.2*m);

	axis -> Draw("histo");
	histos[0] -> DrawNormalized("histo SAME");
	histos[1] -> DrawNormalized("histo SAME");
	histos[2] -> DrawNormalized("histo SAME");
	histos[3] -> DrawNormalized("histo SAME");
	histos[4] -> DrawNormalized("SAME E1");

	leg -> Draw("SAME");
	CMS_lumi(c, 0, 10);

	c -> SaveAs(title + "Signal.png");
	c -> SaveAs(title + "Signal.pdf");



	TCanvas* c2 = new TCanvas();
	c2 -> cd();

	histos[6] -> Add(histos[5]);
	histos[7] -> Add(histos[6]);


	double max2 = histos[4]->GetBinContent(histos[4]->GetMaximumBin());
	if(histos[7]->GetBinContent(histos[7]->GetMaximumBin()) > max2 )
		max2 = histos[7]->GetBinContent(histos[7]->GetMaximumBin());
	histos[7]-> GetYaxis() -> SetRangeUser(0, 1.2*max2);
	histos[7]-> GetYaxis() -> SetTitleSize(0.05);
	histos[7]-> GetYaxis() -> SetTitleFont(42);
	histos[7]-> GetYaxis() -> SetLabelSize(0.045);
	histos[7]-> GetYaxis() -> SetLabelFont(42);
	histos[7]-> GetXaxis() -> SetLabelSize(0);

	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.35, 1, 0.97);
	pad1->SetBottomMargin(0.035);
	pad1->Draw();
	pad1->cd();


	histos[7] -> Draw("histo");
	histos[6] -> Draw("histo SAME");
	histos[5] -> Draw("histo SAME");
	histos[4] -> Draw("SAME E1");
	leg2 -> Draw("SAME");

	c2->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.10, 1, 0.35);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0);
	pad2->Draw();
	pad2->cd();

	TH1F *h = (TH1F*)histos[4]->Clone("h");
	h->SetLineColor(kBlack);
	h->SetMinimum(0.5);  // Define Y ..
	h->SetMaximum(1.5); // .. range
	h->Sumw2();
	h->SetStats(0);      // No statistics on lower plot
	h->Divide(histos[7]);
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
	lumi_sqrtS = "36.73 fb^{-1} (13 TeV)";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

	setTDRStyle();
	gStyle -> SetOptFit(0);
	gStyle -> SetOptStat(0);

	TString HadrFoldMC = "/afs/cern.ch/work/a/abeschi/ttHHadr_v6/";
	TString LeptFoldMC = "/afs/cern.ch/work/a/abeschi/ttHLept_v6/";
	TString HadrFoldData = "/afs/cern.ch/work/a/abeschi/ttHHadr_v6/";
	TString LeptFoldData = "/afs/cern.ch/work/a/abeschi/ttHLept_v6/";

	TChain* ttH;
	TChain* ggH;
	TChain* vbf;
	TChain* vh;
	TChain* data;
	TChain* diphotons;
	TChain* gjet;
	TChain* qcd;

	float lumiFactor = 36.46;
	bool isLept = (argc!=1 ? 0 : 1);

	if(isLept)
	{	cout << "Processing leptonic tag" << endl;
		ttH = new TChain("TTHLeptonicDumper/trees/tree");
		ttH -> Add(LeptFoldMC + "ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_*.root");		
		ggH = new TChain("TTHLeptonicDumper/trees/tree");
		ggH -> Add(LeptFoldMC + "GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8_*.root");		
		vbf = new TChain("TTHLeptonicDumper/trees/tree");
		vbf -> Add(LeptFoldMC + "VBFHToGG_M125_13TeV_amcatnloFXFX_pythia8_*.root");		
		vh = new TChain("TTHLeptonicDumper/trees/tree");
		vh -> Add(LeptFoldMC + "VHToGG_M125_13TeV_amcatnloFXFX_pythia8_*.root");		
		data = new TChain("TTHLeptonicDumper/trees/tree");
		data -> Add(LeptFoldData + "output_DoubleEG_Run2016*.root");		
		diphotons = new TChain("TTHLeptonicDumper/trees/tree");
		diphotons -> Add(LeptFoldMC + "output_DiPhotonJetsBox_MGG*.root");
		gjet = new TChain("TTHLeptonicDumper/trees/tree");
		gjet -> Add(LeptFoldMC + "output_GJet_Pt*.root");		
		qcd = new TChain("TTHLeptonicDumper/trees/tree");
		qcd -> Add(LeptFoldMC + "output_QCD_Pt*.root");

/*		ttH -> Print();
		ggH -> Print();
		vbf -> Print();
		vh -> Print();
		data -> Print();
		diphotons -> Print();
		gjet -> Print();
		qcd -> Print();

		ttH -> SetBranchStatus("*", 1);
		ggH -> SetBranchStatus("*", 1);
		vbf -> SetBranchStatus("*", 1);
		vh -> SetBranchStatus("*", 1);
		data -> SetBranchStatus("*", 1);
		diphotons -> SetBranchStatus("*", 1);
		gjet -> SetBranchStatus("*", 1);
		qcd -> SetBranchStatus("*", 1);
*/
	}

	else
	{	cout << "Processing hadronic tag" << endl;
		ttH = new TChain("TTHHadronicDumper/trees/tree");
		ttH -> Add(HadrFoldMC + "ttHJetToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_*.root");		
		ggH = new TChain("TTHHadronicDumper/trees/tree");
		ggH -> Add(HadrFoldMC + "GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8_*.root");		
		vbf = new TChain("TTHHadronicDumper/trees/tree");
		vbf -> Add(HadrFoldMC + "VBFHToGG_M125_13TeV_amcatnloFXFX_pythia8_*.root");
		vh = new TChain("TTHHadronicDumper/trees/tree");		
		vh -> Add(HadrFoldMC + "VHToGG_M125_13TeV_amcatnloFXFX_pythia8*.root");		
		data = new TChain("TTHHadronicDumper/trees/tree");
		data -> Add(HadrFoldData + "output_DoubleEG_Run2016*.root");		
		diphotons = new TChain("TTHHadronicDumper/trees/tree");
		diphotons -> Add(HadrFoldMC + "output_DiPhotonJetsBox_MGG*.root");
		gjet = new TChain("TTHHadronicDumper/trees/tree");
		gjet -> Add(HadrFoldMC + "output_GJet_Pt*.root");		
		qcd = new TChain("TTHHadronicDumper/trees/tree");
		qcd -> Add(HadrFoldMC + "output_QCD_Pt*.root");	

/*		ttH -> Print();
		ggH -> Print();
		vbf -> Print();
		vh -> Print();
		data -> Print();
		diphotons -> Print();
		gjet -> Print();
		qcd -> Print();
*/


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

	float jet_pt[6] = {0, 0, 0, 0, 0, 0};
	float jet_eta[6] = {0, 0, 0, 0, 0, 0};
	float jet_phi[6] = {0, 0, 0, 0, 0, 0};
	float jet_bdiscriminant[6] = {0, 0, 0, 0, 0, 0};

	float ele_pt[2] = {0, 0};
	float ele_eta[2] = {0, 0};
	float ele_phi[2] = {0, 0};

	float mu_pt[2] = {0, 0};
	float mu_eta[2] = {0, 0};
	float mu_phi[2] = {0, 0};

	float njet = 0;
	float nbjet = 0;


	float Counter[8] = {0., 0., 0., 0., 0., 0., 0., 0.};
	float Entries[8] = {0., 0., 0., 0., 0., 0., 0., 0.};
	TString process[8] = {"ttH", "ggH", "VBF", "VH", "Data", "Diphoton", "Gamma + jets", "QCD" };

	TH1F* nvtx_histo[8];

	TH1F* dipho_sumpt_histo[8];
	TH1F* dipho_mass_histo[8];
	TH1F* dipho_leadPt_histo[8];
	TH1F* dipho_leadEta_histo[8];
	TH1F* dipho_leadPhi_histo[8];
	TH1F* dipho_leadR9_histo[8];
	TH1F* dipho_subleadPt_histo[8];
	TH1F* dipho_subleadEta_histo[8];
	TH1F* dipho_subleadPhi_histo[8];
	TH1F* dipho_subleadR9_histo[8];
	TH1F* dipho_leadIDMVA_histo[8];
	TH1F* dipho_subleadIDMVA_histo[8];
	TH1F* dipho_mva_histo[8];
	TH1F* dipho_leadPtOverM_histo[8];
	TH1F* dipho_subleadPtOverM_histo[8];
	TH1F* jet_pt_histo[8];
	TH1F* jet_eta_histo[8];
	TH1F* jet_phi_histo[8];
	TH1F* jet_bdiscriminant_histo[8];
	TH1F* leadingjet_pt_histo[8];
	TH1F* leadingjet_eta_histo[8];
	TH1F* leadingjet_phi_histo[8];
	TH1F* leadingjet_bdiscriminant_histo[8];
	TH1F* njet_histo[8];
	TH1F* nbjet_histo[8];
	TH1F* nleptons_histo[8];
	TH1F* ele_pt_histo[8];
	TH1F* ele_eta_histo[8];
	TH1F* ele_phi_histo[8];
	TH1F* mu_pt_histo[8];
	TH1F* mu_eta_histo[8];
	TH1F* mu_phi_histo[8];

	for(int i=0; i<8; i++)
	{	
		nvtx_histo[i] = new TH1F(("nvtx_histo"+ std::to_string(i)).c_str(), "; n_{vtx}; Counts", 60, -0.5, 59.5 );
		
		dipho_sumpt_histo[i] = new TH1F(("dipho_sumpt_histo"+ std::to_string(i)).c_str(), "; P_{T}^{#gamma#gamma} (GeV); Counts", 100, 40, 240 );
		dipho_mass_histo[i] = new TH1F(("dipho_mass_histo"+ std::to_string(i)).c_str(), "; m_{#gamma#gamma} (GeV); Counts", 80, 100, 180 );
		dipho_leadPt_histo[i] = new TH1F(("dipho_leadPt_histo"+ std::to_string(i)).c_str(), "; P_{T}^{leading #gamma} (GeV); Counts", 100, 0, 200 );
		dipho_leadEta_histo[i] = new TH1F(("dipho_leadEta_histo"+ std::to_string(i)).c_str(), "; #eta; Counts", 50, -2.5, 2.5 );
		dipho_leadPhi_histo[i] = new TH1F(("dipho_leadPhi_histo"+ std::to_string(i)).c_str(), "; #varphi; Counts", 50, -3.15, 3.15 );
		dipho_leadR9_histo[i] = new TH1F(("dipho_leadR9_histo"+ std::to_string(i)).c_str(), "; R_{9}; Counts", 50, 0.7, 1 );
		dipho_subleadPt_histo[i] = new TH1F(("dipho_subleadPt_histo"+ std::to_string(i)).c_str(), "; P_{T}^{subleading #gamma} (GeV); Counts", 75, 0, 150 );
		dipho_subleadEta_histo[i] = new TH1F(("dipho_subleadEta_histo"+ std::to_string(i)).c_str(), "; #eta; Counts", 50, -2.5, 2.5 );
		dipho_subleadPhi_histo[i] = new TH1F(("dipho_subleadPhi_histo"+ std::to_string(i)).c_str(), "; #varphi; Counts", 50, -3.15, 3.15 );
		dipho_subleadR9_histo[i] = new TH1F(("dipho_subleadR9_histo"+ std::to_string(i)).c_str(), "; R_{9}; Counts", 50, 0.7, 1 );
		dipho_leadIDMVA_histo[i] = new TH1F(("dipho_leadIDMVA_histo"+ std::to_string(i)).c_str(), "; Photon IDMVA; Counts", 50, -1, 1 );
		dipho_subleadIDMVA_histo[i] = new TH1F(("dipho_subleadIDMVA_histo"+ std::to_string(i)).c_str(), "; Photon IDMVA; Counts", 50, -1, 1 );
		dipho_mva_histo[i] = new TH1F(("dipho_mva_histo"+ std::to_string(i)).c_str(), "; Diphoton MVA; Counts", 50, -1, 1 );
		dipho_leadPtOverM_histo[i] = new TH1F(("dipho_leadPtOverM_histo"+ std::to_string(i)).c_str(), "; p_{T}^{lead #gamma}/m_{#gamma #gamma}; Counts", 50, 0.25, 2);
		dipho_subleadPtOverM_histo[i] = new TH1F(("dipho_subleadPtOverM_histo"+ std::to_string(i)).c_str(), "; p_{T}^{lead #gamma}/m_{#gamma #gamma}; Counts", 50, 0.25, 2);
		jet_pt_histo[i] = new TH1F(("jet_pt_histo"+ std::to_string(i)).c_str(), "; Jet p_{T} (GeV); Counts", 100, 0, 200 );
		jet_eta_histo[i] = new TH1F(("jet_eta_histo"+ std::to_string(i)).c_str(), "; Jet #eta; Counts", 50, -2.5, 2.5 );
		jet_phi_histo[i] = new TH1F(("jet_phi_histo"+ std::to_string(i)).c_str(), "; Jet #varphi; Counts", 50, -3.15, 3.15 );
		jet_bdiscriminant_histo[i] = new TH1F(("jet_bdiscriminant_histo"+ std::to_string(i)).c_str(), "; Jet bdiscriminant; Counts", 50, 0, 1 );
		leadingjet_pt_histo[i] = new TH1F(("leadingjet_pt_histo"+ std::to_string(i)).c_str(), "; Jet p_{T} (GeV); Counts", 100, 0, 200 );
		leadingjet_eta_histo[i] = new TH1F(("leadingjet_eta_histo"+ std::to_string(i)).c_str(), "; Jet #eta; Counts", 50, -2.5, 2.5 );
		leadingjet_phi_histo[i] = new TH1F(("leadingjet_phi_histo"+ std::to_string(i)).c_str(), "; Jet #varphi; Counts", 50, -3.15, 3.15 );
		leadingjet_bdiscriminant_histo[i] = new TH1F(("leadingjet_bdiscriminant_histo"+ std::to_string(i)).c_str(), "; Jet bdiscriminant; Counts", 50, 0, 1 );

		njet_histo[i] = new TH1F(("njet_histo"+ std::to_string(i)).c_str(), "; Number of jets; Counts", 16, -0.5, 15.5);
		nbjet_histo[i] = new TH1F(("nbjet_histo"+ std::to_string(i)).c_str(), "; Number of b-jets; Counts", 6, -0.5, 5.5);
		nleptons_histo[i] = new TH1F(("nleptons_histo"+ std::to_string(i)).c_str(), "; Number of leptons; Counts", 6, -0.5, 5.5);

		ele_pt_histo[i] = new TH1F(("ele_pt_histo"+ std::to_string(i)).c_str(), "; Ele p_{T} (GeV); Counts", 100, 0, 200 );
		ele_eta_histo[i] = new TH1F(("ele_eta_histo"+ std::to_string(i)).c_str(), "; Ele #eta; Counts", 50, -2.5, 2.5 );
		ele_phi_histo[i] = new TH1F(("ele_phi_histo"+ std::to_string(i)).c_str(), "; Ele #varphi; Counts", 50, -3.15, 3.15 );
		mu_pt_histo[i] = new TH1F(("mu_pt_histo"+ std::to_string(i)).c_str(), "; Mu p_{T} (GeV); Counts", 100, 0, 200 );
		mu_eta_histo[i] = new TH1F(("mu_eta_histo"+ std::to_string(i)).c_str(), "; Mu #eta; Counts", 50, -2.5, 2.5 );
		mu_phi_histo[i] = new TH1F(("mu_phi_histo"+  std::to_string(i)).c_str(), "; Mu #varphi; Counts", 50, -3.15, 3.15 );
	}

	for(int n=0; n<8; n++)
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
				nentries = data -> GetEntries();
				serviceTree = (TChain*)data -> Clone();
				break;


			case(5):
				nentries = diphotons -> GetEntries();
				serviceTree = (TChain*)diphotons -> Clone();
				break;

			case(6):
				nentries = gjet -> GetEntries();
				serviceTree = (TChain*)gjet -> Clone();
				break;


			case(7):
				nentries = qcd -> GetEntries();
				serviceTree = (TChain*)qcd -> Clone();
				break;


			default:
				nentries = 0;
		}

		cout << "Process: " << process[n] << ", number of event in the tree: " << nentries << endl;

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
			for(int i=1; i<5; i++)
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

		}


		else
		{
			for(int i=1; i<7; i++)
			{	serviceTree -> SetBranchAddress(("jet_pt"+ std::to_string(i)).c_str(), &jet_pt[i-1]);
				serviceTree -> SetBranchAddress(("jet_eta"+ std::to_string(i)).c_str(), &jet_eta[i-1]);
				serviceTree -> SetBranchAddress(("jet_phi"+ std::to_string(i)).c_str(), &jet_phi[i-1]);
				serviceTree -> SetBranchAddress(("jet_bdiscriminant"+ std::to_string(i)).c_str(), &jet_bdiscriminant[i-1]);
			}
			serviceTree -> SetBranchAddress("njet", &njet);
			serviceTree -> SetBranchAddress("nbjet", &nbjet);

		}

		Entries[n] = nentries;

		for(int i=0; i<nentries; i++)
		{	
			serviceTree -> GetEntry(i);
			int nleptons = 0;

			if(dipho_mass<100 || dipho_mass>180) continue;
			if(n==4 && (dipho_mass>115 && dipho_mass<135)) continue;
			if(n==4)
				weight = 1./lumiFactor;

			dipho_mass_histo[n] -> Fill(dipho_mass, weight*lumiFactor);

			if(n>4 && (dipho_mass>115 && dipho_mass<135)) continue;

			Counter[n] += weight*lumiFactor;

			nvtx_histo[n] -> Fill(nvtx, weight*lumiFactor);
			dipho_sumpt_histo[n] -> Fill(dipho_sumpt, weight*lumiFactor);
			dipho_leadPt_histo[n] -> Fill(dipho_leadPt, weight*lumiFactor);
			dipho_leadEta_histo[n] -> Fill(dipho_leadEta, weight*lumiFactor);
			dipho_leadPhi_histo[n] -> Fill(dipho_leadPhi, weight*lumiFactor);
			dipho_leadR9_histo[n] -> Fill(dipho_leadR9, weight*lumiFactor);
			dipho_subleadPt_histo[n] -> Fill(dipho_subleadPt, weight*lumiFactor);
			dipho_subleadEta_histo[n] -> Fill(dipho_subleadEta, weight*lumiFactor);
			dipho_subleadPhi_histo[n] -> Fill(dipho_subleadPhi, weight*lumiFactor);
			dipho_subleadR9_histo[n] -> Fill(dipho_subleadR9, weight*lumiFactor);
			dipho_leadIDMVA_histo[n] -> Fill(dipho_leadIDMVA, weight*lumiFactor);
			dipho_subleadIDMVA_histo[n] -> Fill(dipho_subleadIDMVA, weight*lumiFactor);
			dipho_mva_histo[n] -> Fill(dipho_mva, weight*lumiFactor);
			dipho_leadPtOverM_histo[n] -> Fill(dipho_leadPt/dipho_mass, weight*lumiFactor);
			dipho_subleadPtOverM_histo[n] -> Fill(dipho_subleadPt/dipho_mass, weight*lumiFactor);

			if(jet_pt[0]>0)
			{
				jet_pt_histo[n] -> Fill(jet_pt[0], weight*lumiFactor);
				jet_eta_histo[n] -> Fill(jet_eta[0], weight*lumiFactor);
				jet_phi_histo[n] -> Fill(jet_phi[0], weight*lumiFactor);
				jet_bdiscriminant_histo[n] -> Fill(jet_bdiscriminant[0], weight*lumiFactor);
				leadingjet_pt_histo[n] -> Fill(jet_pt[0], weight*lumiFactor);
				leadingjet_eta_histo[n] -> Fill(jet_eta[0], weight*lumiFactor);
				leadingjet_phi_histo[n] -> Fill(jet_phi[0], weight*lumiFactor);
				leadingjet_bdiscriminant_histo[n] -> Fill(jet_bdiscriminant[0], weight*lumiFactor);
			}

			if(jet_pt[1]>0)
			{
				jet_pt_histo[n] -> Fill(jet_pt[1], weight*lumiFactor);
				jet_eta_histo[n] -> Fill(jet_eta[1], weight*lumiFactor);
				jet_phi_histo[n] -> Fill(jet_phi[1], weight*lumiFactor);
				jet_bdiscriminant_histo[n] -> Fill(jet_bdiscriminant[1], weight*lumiFactor);
			}

			if(jet_pt[2]>0)
			{
				jet_pt_histo[n] -> Fill(jet_pt[2], weight*lumiFactor);
				jet_eta_histo[n] -> Fill(jet_eta[2], weight*lumiFactor);
				jet_phi_histo[n] -> Fill(jet_phi[2], weight*lumiFactor);
				jet_bdiscriminant_histo[n] -> Fill(jet_bdiscriminant[2], weight*lumiFactor);
			}

			if(jet_pt[3]>0)
			{
				jet_pt_histo[n] -> Fill(jet_pt[3], weight*lumiFactor);
				jet_eta_histo[n] -> Fill(jet_eta[3], weight*lumiFactor);
				jet_phi_histo[n] -> Fill(jet_phi[3], weight*lumiFactor);
				jet_bdiscriminant_histo[n] -> Fill(jet_bdiscriminant[3], weight*lumiFactor);
			}


			if(isLept)
			{
				if(ele_pt[0]>0)
				{
					ele_pt_histo[n] -> Fill(ele_pt[0], weight*lumiFactor);
					ele_eta_histo[n] -> Fill(ele_eta[0], weight*lumiFactor);
					ele_phi_histo[n] -> Fill(ele_phi[0], weight*lumiFactor);

					nleptons++;
				}

				if(ele_pt[1]>0)
				{
					ele_pt_histo[n] -> Fill(ele_pt[1], weight*lumiFactor);
					ele_eta_histo[n] -> Fill(ele_eta[1], weight*lumiFactor);
					ele_phi_histo[n] -> Fill(ele_phi[1], weight*lumiFactor);

					nleptons++;
				}

				if(mu_pt[0]>0)
				{
					mu_pt_histo[n] -> Fill(mu_pt[0], weight*lumiFactor);
					mu_eta_histo[n] -> Fill(mu_eta[0], weight*lumiFactor);
					mu_phi_histo[n] -> Fill(mu_phi[0], weight*lumiFactor);

					nleptons++;
				}

				if(mu_pt[1]>0)
				{
					mu_pt_histo[n] -> Fill(mu_pt[1], weight*lumiFactor);
					mu_eta_histo[n] -> Fill(mu_eta[1], weight*lumiFactor);
					mu_phi_histo[n] -> Fill(mu_phi[1], weight*lumiFactor);

					nleptons++;
				}

			}

			else
			{
				if(jet_pt[4]>0)
				{
					jet_pt_histo[n] -> Fill(jet_pt[4], weight*lumiFactor);
					jet_eta_histo[n] -> Fill(jet_eta[4], weight*lumiFactor);
					jet_phi_histo[n] -> Fill(jet_phi[4], weight*lumiFactor);
					jet_bdiscriminant_histo[n] -> Fill(jet_bdiscriminant[4], weight*lumiFactor);
				}

				if(jet_pt[5]>0)
				{
					jet_pt_histo[n] -> Fill(jet_pt[5], weight*lumiFactor);
					jet_eta_histo[n] -> Fill(jet_eta[5], weight*lumiFactor);
					jet_phi_histo[n] -> Fill(jet_phi[5], weight*lumiFactor);
					jet_bdiscriminant_histo[n] -> Fill(jet_bdiscriminant[5], weight*lumiFactor);
				}
			}


			njet_histo[n] -> Fill(njet, weight*lumiFactor);
			nbjet_histo[n] -> Fill(nbjet, weight*lumiFactor);
			nleptons_histo[n] -> Fill(nleptons, weight*lumiFactor);
		}

		delete serviceTree;
	}

	MakePlot(nvtx_histo, "Vertex");
	MakePlot(dipho_sumpt_histo, "DiphotonSumpt");
	MakePlot(dipho_mass_histo, "DiphotonMass");
	MakePlot(dipho_leadPt_histo, "LeadingPhotonPt");
	MakePlot(dipho_leadEta_histo, "LeadingPhotonEta");
	MakePlot(dipho_leadPhi_histo, "LeadingPhotonPhi", 2);
	MakePlot(dipho_leadR9_histo, "LeadingPhotonR9", 1);
	MakePlot(dipho_subleadPt_histo, "SubleadingPhotonPt");
	MakePlot(dipho_subleadEta_histo, "SubleadingPhotonEta");
	MakePlot(dipho_subleadPhi_histo, "SubleadingPhotonPhi", 2);
	MakePlot(dipho_subleadR9_histo, "SubleadingPhotonR9", 1);
	MakePlot(dipho_leadIDMVA_histo, "LeadingPhotonIDMVA", 1);
	MakePlot(dipho_subleadIDMVA_histo, "SubleadingPhotonIDMVA");
	MakePlot(dipho_leadPtOverM_histo, "LeadingPhotonPtOverMass");
	MakePlot(dipho_subleadPtOverM_histo, "SubleadingPhotonPtOverMass");
	MakePlot(dipho_mva_histo, "DiphotonMVA");
	MakePlot(jet_pt_histo, "JetPt");
	MakePlot(jet_eta_histo, "JetEta");
	MakePlot(jet_phi_histo, "JetPhi", 2);
	MakePlot(jet_bdiscriminant_histo, "JetBdiscrimiant");
	MakePlot(leadingjet_pt_histo, "LeadingJetPt");
	MakePlot(leadingjet_eta_histo, "LeadingJetEta");
	MakePlot(leadingjet_phi_histo, "LeadingJetPhi", 2);
	MakePlot(leadingjet_bdiscriminant_histo, "LeadingJetBdiscrimiant");

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
	{	MakePlot(njet_histo, "Njet");
		MakePlot(nbjet_histo, "Nbjet");
	}



	if(isLept)
	{	system("mv *.png ~/www/ttH/Leptonic");
		system("mv *.pdf ~/www/ttH/Leptonic");
	}
	else
	{	system("mv *.png ~/www/ttH/Hadronic");
		system("mv *.pdf ~/www/ttH/Hadronic");
	}


	cout << endl << endl << " *******************************************" << endl << endl;
	for(int i=0; i<8; i++)
	{
		cout << "Process: " << process[i] << " , number of events ";
		if(i>=4)
			cout << "in sidebands -> ";
		else
			cout << " -> ";

		cout << Counter[i] << " ( " << Entries[i] << " unweighted)" << endl;
	}

	cout << endl << " Total events in backgorund MC: " << Counter[5] + Counter[6] + Counter[7] << endl;
	cout << endl  << " *******************************************" << endl << endl;

























}
