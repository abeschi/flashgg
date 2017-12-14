/*
g++ -Wall -o EfficiencySemiLeptonic `root-config --cflags --glibs` -L $ROOTSYS/lib -lRooFitCore -lFoam -lMinuit -lMathMore CMS_lumi.C tdrstyle.C EfficiencySemiLeptonic.cpp
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
#include "TLorentzVector.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"

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



void MakePlot(TH1F* histo, TString title)
{
	TCanvas* c = new TCanvas;
	c -> cd();

	histo -> SetLineWidth(2);
	histo -> SetLineColor(kBlue);
	histo -> GetYaxis() -> SetRangeUser(0, 1.2*histo -> GetMaximum());
	histo -> Draw("histo");

	CMS_lumi(c, 0, 0);

	c -> SaveAs(title + ".pdf");
	c -> SaveAs(title + ".png");

	return;
}


void MakeEfficiencyPlot(TEfficiency* eff1, TEfficiency* eff2, float integratedEfficiency1, float integratedEfficiency2, TString title)
{

	eff1 -> SetMarkerStyle(20);
	eff1 -> SetMarkerColor(kRed+1);
	eff1 -> SetLineColor(kRed+1);
	eff1 -> SetFillStyle(0);

	eff2 -> SetMarkerStyle(21);
	eff2 -> SetMarkerColor(kGreen+2);
	eff2 -> SetLineColor(kGreen+1);
	eff2 -> SetFillStyle(0);

	TCanvas* c = new TCanvas;
	c -> cd();

	eff1 -> Draw("");
	eff2 -> Draw("SAME");

	TLegend* leg = new TLegend(0.18, 0.4, 0.5, 0.6);
	leg -> AddEntry(eff1, "Diphoton vertex", "p");
	leg -> AddEntry(eff2, "Vertex 0", "p");
	leg -> Draw("SAME");

	c -> Update();
	eff1 -> GetPaintedGraph() -> GetYaxis() -> SetRangeUser(0.95, 1.01);

	stringstream ss1;
	ss1 << "Diphoton vertex efficiency: " << std::setprecision(4) << integratedEfficiency1*100 << "%" << endl;
	TString string1 = ss1.str();

	TLatex* tex1 = new TLatex(0.18, 0.27, string1);
	tex1 -> SetNDC();
	tex1 -> SetTextColor(kRed+1);
	tex1 -> Draw("SAME");

	stringstream ss2;
	ss2 << "Vertex 0 efficiency: " << std::setprecision(4) << integratedEfficiency2*100 << "%" << endl;
	TString string2 = ss2.str();

	TLatex* tex2 = new TLatex(0.18, 0.2, string2);
	tex2 -> SetNDC();
	tex2 -> SetTextColor(kGreen+2);
	tex2 -> Draw("SAME");

	CMS_lumi(c, 0, 0);

	c -> SaveAs(title + ".pdf");
	c -> SaveAs(title + ".png");

	return;
}


void MakePlot2D(TH2F* histo, TString title)
{
	TCanvas* c = new TCanvas;
	c -> cd();

	histo -> Draw("COLZ");

	CMS_lumi(c, 0, 0);

	c -> SaveAs(title + ".pdf");
	c -> SaveAs(title + ".png");

	return;
}

void Make2Plot(TH1F* histo1, TH1F* histo2, TString title, TString l1, TString l2, bool Normalize=0)
{
	TCanvas* c = new TCanvas;
	c -> cd();

	histo1 -> SetLineWidth(3);
	histo1 -> SetLineColor(kRed+1);

	histo2 -> SetLineWidth(3);
	histo2 -> SetLineColor(kGreen+2);

	float Max = 0;
	if(!Normalize)
		Max = max(histo1->GetMaximum(), histo2->GetMaximum());
	else
		Max = max(histo1->GetMaximum()/histo1->Integral(), histo2->GetMaximum()/histo2->Integral());

	histo1 -> GetYaxis() -> SetRangeUser(0, 1.1*Max);
	if(Normalize)
	{
		histo1 -> DrawNormalized("histo");
		histo2 -> DrawNormalized("histo SAME");	
	}	

	else
	{
		histo1 -> Draw("histo");
		histo2 -> Draw("histo SAME");	
	}	

	TLegend* leg = new TLegend(0.5, 0.82, 0.9, 0.92);
	leg -> AddEntry(histo1, l1, "l");
	leg -> AddEntry(histo2, l2, "l");
	leg -> Draw("SAME");

	CMS_lumi(c, 0, 0);

	c -> SaveAs(title + ".pdf");
	c -> SaveAs(title + ".png");

	return;
}



void MakeAllPlots(const int n, TH1F** histo, TString title, bool Normalize=0)
{
	TCanvas* c = new TCanvas;
	c -> cd();

	int Color[n] = {kOrange-7, kOrange-3, kRed+2, kRed, kOrange, kYellow, kSpring, kGreen+1, kTeal, kAzure+5, kAzure, kViolet+2, kRed+1};

	for(int i=0; i<n; i++)
	{	histo[i] -> SetLineWidth(3);
		histo[i] -> SetLineColor(Color[i]);
		histo[i] -> SetFillStyle(0);
	}


	TCanvas* allMVA = new TCanvas();
	allMVA -> cd();
	
	histo[0] -> Draw("histo");
	for(int i=1; i<n; i++)
		histo[i] -> Draw("histo SAME");

	TLegend* leg1 = new TLegend(0.55, 0.65, 0.9, 0.9);
	for(int i=0; i<n; i++)
		leg1 -> AddEntry(histo[i], histo[i]->GetTitle(), "l");
	leg1 -> Draw("SAME");

	CMS_lumi(allMVA, 0, 0);

	float integratedEfficiency1 = histo[n-1]->GetSumOfWeights()/histo[0]->GetSumOfWeights();
	stringstream ss1;
	ss1 << "#varepsilon: " << std::setprecision(4) << integratedEfficiency1*100 << "%" << endl;
	TString string1 = ss1.str();

	TLatex* tex1 = new TLatex(0.65, 0.50, string1);
	tex1 -> SetNDC();
	tex1 -> Draw("SAME");


	allMVA -> SaveAs(title + "Fully.pdf");
	allMVA -> SaveAs(title + "Fully.png");




	TCanvas* RecoMVA = new TCanvas();
	RecoMVA -> cd();
	
	histo[3] -> Draw("histo");
	for(int i=4; i<n; i++)
		histo[i] -> Draw("histo SAME");

	TLegend* leg3 = new TLegend(0.6, 0.65, 0.9, 0.9);
	for(int i=3; i<n; i++)
		leg3 -> AddEntry(histo[i], histo[i]->GetTitle(), "l");
	leg3 -> Draw("SAME");

	CMS_lumi(RecoMVA, 0, 0);

	RecoMVA -> SaveAs(title + "RecoFully.pdf");
	RecoMVA -> SaveAs(title + "RecoFully.png");



	return;
}








int main(int argc, char *argv[])
{
	writeExtraText = true;       // if extra text
	extraText  = "Preliminary Simulation";  // default extra text is "Preliminary"
	lumi_sqrtS = "(13 TeV)";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

	setTDRStyle();
	gStyle -> SetOptFit(0);
	gStyle -> SetOptStat(0);

	float bJetLooseThreshold = 0.5426;
	float bJetMediumThreshold = 0.8484;
	float jetPtThreshold = 25.;
	float jetEtaThreshold = 2.4;
	float MZ = 91.187;
	bool excludeTaus = false;
	bool onlyTaus = false;

	TString SignalFold = "/afs/cern.ch/work/a/abeschi/ttHEfficiency_v5/";

	TChain* ttH;
	
	float lumiFactor = 35.9;

	ttH = new TChain("TTHFullyLeptonicEfficiencyDumper/trees/tth_13TeV_all");
	ttH -> Add(SignalFold + "ttHJetToGG_M125_13TeV_*.root");		

	int nvtx = 0;
	float weight = 0;
	
	float dipho_sumpt = 0;
	float dipho_mass = 0;
	float dipho_leadPt = 0;
	float dipho_leadEnergy = 0;
	float dipho_leadEta = 0;
	float dipho_leadPhi = 0;
	float dipho_leadR9 = 0;
	float dipho_lead_hoe = 0;
	float dipho_lead_ptoM = 0;
	float dipho_leadEgChargedHadronIso = 0;
	float dipho_subleadPt = 0;
	float dipho_subleadEnergy = 0;
	float dipho_subleadEta = 0;
	float dipho_subleadPhi = 0;
	float dipho_subleadR9 = 0;
	float dipho_sublead_hoe = 0;
	float dipho_sublead_ptoM = 0;
	float dipho_subleadEgChargedHadronIso = 0;
	float dipho_leadIDMVA = 0;
	float dipho_subleadIDMVA = 0;
	float dipho_mva = 0;
	float dipho_vertexZ = 0;
	float dipho_vertexSigmaZ = 0;

	float passHLT = 0;
	float passPreselection = 0;
	
	float jet_pt[5] = {0., 0., 0., 0., 0.};
	float jet_eta[5] = {0., 0., 0., 0., 0.};
	float jet_phi[5] = {0., 0., 0., 0., 0.};
	float jet_bdiscriminant[5] = {0., 0., 0., 0., 0.};
	float jet_energy[5] = {0., 0., 0., 0., 0.};

	float lepton_pt[4] = {0., 0., 0., 0.};
	float lepton_eta[4] = {0., 0., 0., 0.};
	float lepton_phi[4] = {0., 0., 0., 0.};

	float MetPt = 0.;
	float MetPhi = 0.;
	float njet = 0;
	float nbjetLoose = 0;
	float nbjetMedium = 0;
	float nbjetTight = 0;

	float GenLeadPhoton_pt = 0.;
	float GenLeadPhoton_eta = 0.;
	float GenLeadPhoton_phi = 0.;
	float GenLeadPhoton_energy = 0.;
	float GenSubleadPhoton_pt = 0.;
	float GenSubleadPhoton_eta = 0.;
	float GenSubleadPhoton_phi = 0.;
	float GenSubleadPhoton_energy = 0.;
	float HiggsVtxZ = 0.;
	float diphotonSize = 0.;

	float GenJet_pt[4] = {0., 0., 0., 0.};
	float GenJet_eta[4] = {0., 0., 0., 0.};
	float GenJet_phi[4] = {0., 0., 0., 0.};
	float GenJet_energy[4] = {0., 0., 0., 0.};
	float GenLepton_pt[2] = {0., 0.};
	float GenLepton_eta[2] = {0., 0.};
	float GenLepton_phi[2] = {0., 0.};
	float GenLepton_energy[2] = {0., 0.};
	float GenLepton_pdgId[2] = {0., 0.};
	float NGenJets = 0;
	float NGenJetsEta2p4 = 0;
	float Vtx0Z = 0;
	float Vtx0ZSigma = 0;

	ttH -> SetBranchAddress("nvtx", &nvtx);
	ttH -> SetBranchAddress("weight", &weight);
	ttH -> SetBranchAddress("dipho_sumpt", &dipho_sumpt);
	ttH -> SetBranchAddress("dipho_mass", &dipho_mass);
	ttH -> SetBranchAddress("dipho_leadPt", &dipho_leadPt);
	ttH -> SetBranchAddress("dipho_leadEnergy", &dipho_leadEnergy);
	ttH -> SetBranchAddress("dipho_leadEta", &dipho_leadEta);
	ttH -> SetBranchAddress("dipho_leadPhi", &dipho_leadPhi);
	ttH -> SetBranchAddress("dipho_leadR9", &dipho_leadR9);
	ttH -> SetBranchAddress("dipho_lead_hoe", &dipho_lead_hoe);
	ttH -> SetBranchAddress("dipho_lead_ptoM", &dipho_lead_ptoM);
	ttH -> SetBranchAddress("dipho_leadEgChargedHadronIso", &dipho_leadEgChargedHadronIso);
	ttH -> SetBranchAddress("dipho_subleadPt", &dipho_subleadPt);
	ttH -> SetBranchAddress("dipho_subleadEnergy", &dipho_subleadEnergy);
	ttH -> SetBranchAddress("dipho_subleadEta", &dipho_subleadEta);
	ttH -> SetBranchAddress("dipho_subleadPhi", &dipho_subleadPhi);
	ttH -> SetBranchAddress("dipho_subleadR9", &dipho_subleadR9);
	ttH -> SetBranchAddress("dipho_sublead_hoe", &dipho_sublead_hoe);
	ttH -> SetBranchAddress("dipho_sublead_ptoM", &dipho_sublead_ptoM);
	ttH -> SetBranchAddress("dipho_subleadEgChargedHadronIso", &dipho_subleadEgChargedHadronIso);
	ttH -> SetBranchAddress("dipho_leadIDMVA", &dipho_leadIDMVA);
	ttH -> SetBranchAddress("dipho_subleadIDMVA", &dipho_subleadIDMVA);
	ttH -> SetBranchAddress("dipho_mva", &dipho_mva);
	ttH -> SetBranchAddress("dipho_vertexZ", &dipho_vertexZ);
	ttH -> SetBranchAddress("dipho_vertexSigmaZ", &dipho_vertexSigmaZ);

	ttH -> SetBranchAddress("passHLT", &passHLT);
	ttH -> SetBranchAddress("passPreselection", &passPreselection);

	ttH -> SetBranchAddress("GenLeadPhton_pt", &GenLeadPhoton_pt);
	ttH -> SetBranchAddress("GenLeadPhton_eta", &GenLeadPhoton_eta);
	ttH -> SetBranchAddress("GenLeadPhton_phi", &GenLeadPhoton_phi);
	ttH -> SetBranchAddress("GenLeadPhton_energy", &GenLeadPhoton_energy);
	ttH -> SetBranchAddress("GenSubleadPhton_pt", &GenSubleadPhoton_pt);
	ttH -> SetBranchAddress("GenSubleadPhton_eta", &GenSubleadPhoton_eta);
	ttH -> SetBranchAddress("GenSubleadPhton_phi", &GenSubleadPhoton_phi);
	ttH -> SetBranchAddress("GenSubleadPhton_energy", &GenSubleadPhoton_energy);
	ttH -> SetBranchAddress("HiggsVtxZ", &HiggsVtxZ);
	ttH -> SetBranchAddress("diphotonSize", &diphotonSize);

	for(int i=1; i<6; i++)
	{	ttH -> SetBranchAddress(("jet_pt"+ std::to_string(i)).c_str(), &jet_pt[i-1]);
		ttH -> SetBranchAddress(("jet_eta"+ std::to_string(i)).c_str(), &jet_eta[i-1]);
		ttH -> SetBranchAddress(("jet_phi"+ std::to_string(i)).c_str(), &jet_phi[i-1]);
		ttH -> SetBranchAddress(("jet_bdiscriminant"+ std::to_string(i)).c_str(), &jet_bdiscriminant[i-1]);
		ttH -> SetBranchAddress(("jet_energy"+ std::to_string(i)).c_str(), &jet_energy[i-1]);
	}

	for(int i=1; i<3; i++)
	{	ttH -> SetBranchAddress(("mu_pt"+ std::to_string(i)).c_str(), &lepton_pt[i-1]);
		ttH -> SetBranchAddress(("mu_eta"+ std::to_string(i)).c_str(), &lepton_eta[i-1]);
		ttH -> SetBranchAddress(("mu_phi"+ std::to_string(i)).c_str(), &lepton_phi[i-1]);
		ttH -> SetBranchAddress(("ele_pt"+ std::to_string(i)).c_str(), &lepton_pt[i+1]);
		ttH -> SetBranchAddress(("ele_eta"+ std::to_string(i)).c_str(), &lepton_eta[i+1]);
		ttH -> SetBranchAddress(("ele_phi"+ std::to_string(i)).c_str(), &lepton_phi[i+1]);
	}

	ttH -> SetBranchAddress("MetPt", &MetPt);
	ttH -> SetBranchAddress("MetPhi", &MetPhi);
	ttH -> SetBranchAddress("njet", &njet);
	ttH -> SetBranchAddress("nbjetLoose", &nbjetLoose);
	ttH -> SetBranchAddress("nbjetMedium", &nbjetMedium);
	ttH -> SetBranchAddress("nbjetTight", &nbjetTight);
	ttH -> SetBranchAddress("NGenJets", &NGenJets);
	ttH -> SetBranchAddress("NGenJetsEta2p4", &NGenJetsEta2p4);
	ttH -> SetBranchAddress("Vtx0Z", &Vtx0Z);
	ttH -> SetBranchAddress("Vtx0ZSigma", &Vtx0ZSigma);

	for(int i=1; i<4; i++)
	{
		ttH -> SetBranchAddress(("GenJet_pt"+ std::to_string(i)).c_str(), &GenJet_pt[i-1]);
		ttH -> SetBranchAddress(("GenJet_eta"+ std::to_string(i)).c_str(), &GenJet_eta[i-1]);
		ttH -> SetBranchAddress(("GenJet_phi"+ std::to_string(i)).c_str(), &GenJet_phi[i-1]);
		ttH -> SetBranchAddress(("GenJet_energy"+ std::to_string(i)).c_str(), &GenJet_energy[i-1]);
	}

	for(int i=1; i<2; i++)
	{
		ttH -> SetBranchAddress(("GenLepton_pt"+ std::to_string(i)).c_str(), &GenLepton_pt[i-1]);
		ttH -> SetBranchAddress(("GenLepton_eta"+ std::to_string(i)).c_str(), &GenLepton_eta[i-1]);
		ttH -> SetBranchAddress(("GenLepton_phi"+ std::to_string(i)).c_str(), &GenLepton_phi[i-1]);
		ttH -> SetBranchAddress(("GenLepton_energy"+ std::to_string(i)).c_str(), &GenLepton_energy[i-1]);
		ttH -> SetBranchAddress(("GenLepton_pdgId"+ std::to_string(i)).c_str(), &GenLepton_pdgId[i-1]);
	}

	int nentries = ttH -> GetEntries();
//	int nentries = 100;

/*	
	Gen
	Gen + geometric acceptance
	Reco -> acceptance*reco efficiency
	+ HLT
	+ preselections
	+ veto on reco lepton
	  - MVA (3 jets, bLoose, tthMVA)
	  - cut based (5 jets, bMedium)
	+ Diphoton MVA
*/

	const int n = 13;

	float efficiencyDiphotonVtx = 0.;
	float efficiencyVtx0 = 0.;
	float totalEvents = 0.;

	TH1F* LeadingEta = new TH1F("LeadingEta", "; #eta^{lead #gamma}; Counts", 100 , -5, 5);
	TH1F* LeadingEtaCut = new TH1F("LeadingEtaCut", "; #eta^{lead #gamma}; Counts", 100 , -5, 5);
	TH1F* SubleadingEta = new TH1F("SubleadingEta", "; #eta^{sublead #gamma}; Counts", 100 , -5, 5);
	TH1F* SubleadingEtaCut = new TH1F("SubleadingEtaCut", "; #eta^{sublead #gamma}; Counts", 100 , -5, 5);

	TH2F* EtaPhotonDistribution = new TH2F("EtaPhotonDistribution", "; #eta^{lead #gamma}; #eta^{sublead #gamma}; Counts", 100 , -5, 5., 100 , -5, 5);

	double PtBins[15] = {0.,100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 230., 250., 350.};
	TEfficiency* VtxVsPt = new TEfficiency("VtxVsPt", "; P_{T}^{#gamma #gamma} (GeV); #varepsilon", 14, PtBins);
	TEfficiency* VtxVsPt0 = new TEfficiency("VtxVsPt0", "; P_{T}^{#gamma #gamma} (GeV); #varepsilon", 14, PtBins);

	double VtxBins[12] = {0.,5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 60.,};
	TEfficiency* VtxVsPileUp = new TEfficiency("VtxVsPileUp", "; n_{vtx}; Counts", 11, VtxBins);
	TEfficiency* VtxVsPileUp0 = new TEfficiency("VtxVsPileUp0", "; n^{vtx}; Counts", 11, VtxBins);

	TH1F* DistanceFromTrueVtx = new TH1F("DistanceFromTrueVtx", "; (z^{true} - z^{reco})/#sigma; couts", 100, 0, 5);
	TH1F* DistanceFromTrueVtx0 = new TH1F("DistanceFromTrueVtx0", "; (z^{true} - z^{reco})/#sigma; couts", 100, 0, 5);

	TH1F* mass[n];
	TH1F* eta[n];
	TH1F* pt[n];
	TString title[n] = {"Generated photons", "Photons not in EB-EE crack", "Photons within acceptance", "Reconstructed photons", "pass HLT", "pass Preselection", "p_{T} photons > 33 (25) GeV", "p_{T} photons > 0.3(0.25)*m_{#gamma #gamma}","p_{T} photons > 0.5(0.25)*m_{#gamma #gamma}", "At least one lepton", "At least 2 jets", "At least 1 Medium b-jet", "Diphoton MVA > -0.4"};

	for(int i=0; i<n; i++)
	{
		mass[i] = new TH1F(("mass" + to_string(i)).c_str(), "; m_{#gamma #gamma} (GeV); Counts", 80 , 115, 135);
		eta[i] = new TH1F(("eta" + to_string(i)).c_str(), "; #eta^{#gamma #gamma}; Counts", 100 , -5., 5.);
		pt[i] = new TH1F(("pt" + to_string(i)).c_str(), "; p_{T}^{#gamma #gamma} (GeV); Counts", 300 , 50., 350.);

		mass[i] -> Sumw2();
		eta[i] -> Sumw2();
		pt[i] -> Sumw2();

		mass[i] -> SetTitle(title[i]);
		eta[i] -> SetTitle(title[i]);
		pt[i] -> SetTitle(title[i]);
	}



	for(int ev=0; ev<nentries; ev++)
	{	
		ttH -> GetEntry(ev);
		if(ev%10000==0) cout << "Processing " << ev << "th event out of " << nentries << "\r" << flush;

		if(excludeTaus && (abs(GenLepton_pdgId[0]==15) || abs(GenLepton_pdgId[1]==15))) continue;
		if(onlyTaus && !(abs(GenLepton_pdgId[0]==15) || abs(GenLepton_pdgId[1]==15))) continue;


		//Invatiant mass of genPhotons
		TLorentzVector* Gen1 = new TLorentzVector();
		Gen1 -> SetPtEtaPhiE(GenLeadPhoton_pt, GenLeadPhoton_eta, GenLeadPhoton_phi, GenLeadPhoton_energy);
		TLorentzVector* Gen2 = new TLorentzVector();
		Gen2 -> SetPtEtaPhiE(GenSubleadPhoton_pt, GenSubleadPhoton_eta, GenSubleadPhoton_phi, GenSubleadPhoton_energy);

		TLorentzVector* Reco1 = new TLorentzVector();
		Reco1 -> SetPtEtaPhiE(dipho_leadPt, dipho_leadEta, dipho_leadPhi, dipho_leadEnergy);
		TLorentzVector* Reco2 = new TLorentzVector();
		Reco2 -> SetPtEtaPhiE(dipho_subleadPt, dipho_subleadEta, dipho_subleadPhi, dipho_subleadEnergy);

		int njets = 0;
		int nbjetsLoose = 0;
		int nbjetsMedium = 0;

		for(int j=0; j<4; j++)
		{	for(int i=0; i<9; i++)
			{
				float dr = compute_R(lepton_eta[j], jet_eta[i], lepton_phi[j], jet_phi[i]);
				if(jet_pt[i]<jetPtThreshold || abs(jet_eta[i])>jetEtaThreshold || dr<0.4) continue;
				njets++;
				if(jet_bdiscriminant[i]>bJetMediumThreshold)
					nbjetsMedium++;
				if(jet_bdiscriminant[i]>bJetLooseThreshold)
					nbjetsLoose++;
			}
		}

		int nLeptons = 0;
		if(lepton_pt[0] > 20.) nLeptons++;
		if(lepton_pt[1] > 20.) nLeptons++;
		float mass1 = sqrt(2*lepton_pt[2]*dipho_leadPt*( ( cosh(lepton_eta[2]-dipho_leadEta) - cos(lepton_phi[2] - dipho_leadPhi) ) ));
		float mass2 = sqrt(2*lepton_pt[2]*dipho_subleadPt*( ( cosh(lepton_eta[2]-dipho_subleadEta) - cos(lepton_phi[2] - dipho_subleadPhi) ) ));
		if(lepton_pt[2]>20 && (mass1-MZ)>5 && (mass2-MZ)>5.) nLeptons++;
		mass1 = sqrt(2*lepton_pt[3]*dipho_leadPt*( ( cosh(lepton_eta[3]-dipho_leadEta) - cos(lepton_phi[3] - dipho_leadPhi) ) ));
		mass2 = sqrt(2*lepton_pt[3]*dipho_subleadPt*( ( cosh(lepton_eta[3]-dipho_subleadEta) - cos(lepton_phi[3] - dipho_subleadPhi) ) ));
		if(lepton_pt[3]>20 && (mass1-MZ)>5 && (mass2-MZ)>5.) nLeptons++;

		LeadingEta -> Fill(GenLeadPhoton_eta, weight*lumiFactor);
		SubleadingEta -> Fill(GenSubleadPhoton_eta, weight*lumiFactor);
		EtaPhotonDistribution -> Fill(GenLeadPhoton_eta, GenSubleadPhoton_eta, weight*lumiFactor);

		mass[0] -> Fill((*Gen1 + *Gen2).M(), weight*lumiFactor);
		eta[0] -> Fill((*Gen1 + *Gen2).Eta(), weight*lumiFactor);
		pt[0] -> Fill(Gen1->Pt() + Gen2->Pt(), weight*lumiFactor);


		if( (abs(GenLeadPhoton_eta)<1.4442 || abs(GenLeadPhoton_eta)>1.556) && (abs(GenSubleadPhoton_eta)<1.4442 || abs(GenSubleadPhoton_eta)>1.556) )
		{	mass[1] -> Fill((*Gen1 + *Gen2).M(), weight*lumiFactor);
			eta[1] -> Fill((*Gen1 + *Gen2).Eta(), weight*lumiFactor);
			pt[1] -> Fill(Gen1->Pt() + Gen2->Pt(), weight*lumiFactor);
		}

		if( (abs(GenLeadPhoton_eta)<1.4442 || abs(GenLeadPhoton_eta)>1.556) && (abs(GenSubleadPhoton_eta)<1.4442 || abs(GenSubleadPhoton_eta)>1.556) && (abs(GenLeadPhoton_eta)<2.4 && abs(GenSubleadPhoton_eta)<2.4) )
		{	mass[2] -> Fill((*Gen1 + *Gen2).M(), weight*lumiFactor);
			eta[2] -> Fill((*Gen1 + *Gen2).Eta(), weight*lumiFactor);
			pt[2] -> Fill(Gen1->Pt() + Gen2->Pt(), weight*lumiFactor);
			LeadingEtaCut -> Fill(GenLeadPhoton_eta, weight*lumiFactor);
			SubleadingEtaCut -> Fill(GenSubleadPhoton_eta, weight*lumiFactor);
		}

		if(dipho_leadPt>0)
		{	mass[3] -> Fill(dipho_mass, weight*lumiFactor);
			eta[3] -> Fill((*Reco1+*Reco2).Eta(), weight*lumiFactor);
			pt[3] -> Fill(dipho_sumpt, weight*lumiFactor);
		}

		if(dipho_leadPt>0 && passHLT==1)
		{	mass[4] -> Fill(dipho_mass, weight*lumiFactor);
			eta[4] -> Fill((*Reco1+*Reco2).Eta(), weight*lumiFactor);
			pt[4] -> Fill(dipho_sumpt, weight*lumiFactor);
		}

		if(dipho_leadPt>0 && passHLT==1 && passPreselection==1)
		{	mass[5] -> Fill(dipho_mass, weight*lumiFactor);
			eta[5] -> Fill((*Reco1+*Reco2).Eta(), weight*lumiFactor);
			pt[5] -> Fill(dipho_sumpt, weight*lumiFactor);
		}

		if(dipho_leadPt>0 && passHLT==1 && passPreselection==1 && dipho_leadPt>100/3. && dipho_subleadPt>25.)
		{	mass[6] -> Fill(dipho_mass, weight*lumiFactor);
			eta[6] -> Fill((*Reco1+*Reco2).Eta(), weight*lumiFactor);
			pt[6] -> Fill(dipho_sumpt, weight*lumiFactor);
		}

		if(dipho_leadPt>0 && passHLT==1 && passPreselection==1 && dipho_leadPt>0.33*dipho_mass && dipho_subleadPt>0.25*dipho_mass)
		{	mass[7] -> Fill(dipho_mass, weight*lumiFactor);
			eta[7] -> Fill((*Reco1+*Reco2).Eta(), weight*lumiFactor);
			pt[7] -> Fill(dipho_sumpt, weight*lumiFactor);
		}

		if(dipho_leadPt>0 && passHLT==1 && passPreselection==1 && dipho_leadPt>0.5*dipho_mass && dipho_subleadPt>0.25*dipho_mass)
		{	mass[8] -> Fill(dipho_mass, weight*lumiFactor);
			eta[8] -> Fill((*Reco1+*Reco2).Eta(), weight*lumiFactor);
			pt[8] -> Fill(dipho_sumpt, weight*lumiFactor);
		}

		if(dipho_leadPt>0 &&passHLT==1 && passPreselection==1 && dipho_leadPt>0.5*dipho_mass && dipho_subleadPt>0.25*dipho_mass && nLeptons>0)
		{	mass[9] -> Fill(dipho_mass, weight*lumiFactor);
			eta[9] -> Fill((*Reco1+*Reco2).Eta(), weight*lumiFactor);
			pt[9] -> Fill(dipho_sumpt, weight*lumiFactor);
		}
				
		if(dipho_leadPt>0 &&passHLT==1 && passPreselection==1 && dipho_leadPt>0.5*dipho_mass && dipho_subleadPt>0.25*dipho_mass && nLeptons>0 && njets>=2)
		{	mass[10] -> Fill(dipho_mass, weight*lumiFactor);
			eta[10] -> Fill((*Reco1+*Reco2).Eta(), weight*lumiFactor);
			pt[10] -> Fill(dipho_sumpt, weight*lumiFactor);
		}

		if(dipho_leadPt>0 &&passHLT==1 && passPreselection==1 && dipho_leadPt>0.5*dipho_mass && dipho_subleadPt>0.25*dipho_mass && nLeptons>0 && njets>=2 && nbjetsMedium>=1)
		{	mass[11] -> Fill(dipho_mass, weight*lumiFactor);
			eta[11] -> Fill((*Reco1+*Reco2).Eta(), weight*lumiFactor);
			pt[11] -> Fill(dipho_sumpt, weight*lumiFactor);
		}

		if(dipho_leadPt>0 &&passHLT==1 && passPreselection==1 && dipho_leadPt>0.5*dipho_mass && dipho_subleadPt>0.25*dipho_mass && nLeptons>0 && njets>=2 && nbjetsMedium>=1 && dipho_mva>-0.4)
		{	mass[12] -> Fill(dipho_mass, weight*lumiFactor);
			eta[12] -> Fill((*Reco1+*Reco2).Eta(), weight*lumiFactor);
			pt[12] -> Fill(dipho_sumpt, weight*lumiFactor);
		}

		if(dipho_leadPt>0)
		{	bool DiphoVtxOk = (HiggsVtxZ - dipho_vertexZ)<1 ? 1: 0;
			bool Vtx0Ok = (HiggsVtxZ - Vtx0Z)<1 ? 1: 0;

			VtxVsPt -> FillWeighted(DiphoVtxOk, weight*lumiFactor, dipho_sumpt);
			VtxVsPileUp -> FillWeighted(DiphoVtxOk, weight*lumiFactor, nvtx);
			VtxVsPt0-> FillWeighted(Vtx0Ok, weight*lumiFactor, dipho_sumpt);
			VtxVsPileUp0 -> FillWeighted(Vtx0Ok, weight*lumiFactor, nvtx);

			DistanceFromTrueVtx -> Fill((HiggsVtxZ - dipho_vertexZ)/dipho_vertexSigmaZ, weight*lumiFactor);
			DistanceFromTrueVtx0 -> Fill((HiggsVtxZ - Vtx0Z)/Vtx0ZSigma, weight*lumiFactor);

			totalEvents += weight*lumiFactor;
			efficiencyDiphotonVtx += weight*lumiFactor*DiphoVtxOk;
			efficiencyVtx0 += weight*lumiFactor*Vtx0Ok;
		}
	}

	cout << "Processing of " << nentries << " events out of " << nentries << " completed" << endl;
	efficiencyDiphotonVtx = efficiencyDiphotonVtx/totalEvents;
	efficiencyVtx0 = efficiencyVtx0/totalEvents;



	float events[n];
	float percentage[n];

	for(int i=0; i<n; i++)
	{
		events[i] = mass[i]->GetSumOfWeights();
		percentage[i] = mass[i]->GetSumOfWeights()/mass[0]->GetSumOfWeights()*100;
	}

	string name[n] = { "Genrated events                          : ", "Not in EB-EE gap                         : ", "Geometric acceptance                     : ", "Acceptance and reconstruction efficiency : ", "HLT efficiency                           : ", "Preselection                             : ", "Photons pT > 33(25) GeV                  : ", "Photons pT > 0.33(0.25) diphoton mass    : ", "Photons pT > 0.5(0.25) diphoton mass     : ", "At least one lepton                      : ", "At least 2 jets                          : ", "At least 1 Medium b-Jet                  : ", "Diphoton MVA > -0.4                      : "};

	cout << endl << endl << "|||--------------------------------------------------------------------------------------|||" << endl;
	cout <<                 "                                        EFFICIENCY                                         " << endl;

	for(int i=0; i<n; i++)
		cout << name[i] << std::setprecision(3) << events[i] << ", " << percentage[i] << "% of the total" << endl;

	cout << endl << "|||--------------------------------------------------------------------------------------|||" << endl;
	cout << endl << endl << endl << endl << endl;





	MakePlot(LeadingEta, "LeadingPhotonEta");
	MakePlot(LeadingEtaCut, "LeadingPhotonEtaCut");
	MakePlot(SubleadingEta, "SubleadingPhotonEta");
	MakePlot(SubleadingEtaCut, "SubleadingPhotonEtaCut");

	Make2Plot(DistanceFromTrueVtx, DistanceFromTrueVtx0, "DistanceFromTrueVtx", "Diphoton Vertex0", "Vertex 0");
	Make2Plot(LeadingEta, LeadingEtaCut, "LeadingPhotonEtaCompare",  "No eta cut", "|#eta|<2.4 and not in EB-EE gap");
	Make2Plot(SubleadingEta, SubleadingEtaCut, "SubleadingPhotonEtaCompare", "No eta cut", "|#eta|<2.4 and not in EB-EE gap");

	MakePlot2D(EtaPhotonDistribution, "EtaPhotonDistribution");

	MakeEfficiencyPlot(VtxVsPt, VtxVsPt0, efficiencyDiphotonVtx, efficiencyVtx0, "VertexEfficiencyVsHiggsPt");
	MakeEfficiencyPlot(VtxVsPileUp, VtxVsPileUp0, efficiencyDiphotonVtx, efficiencyVtx0, "VertexEfficiencyVsPileUp");

	MakeAllPlots(n, mass, "DiphtonMass");
	MakeAllPlots(n, eta, "DiphotonEta");
	MakeAllPlots(n, pt, "DiphotonPt");
	

	system("mv *.png ~/www/ttH/Full2016Dataset/Leptonic/Efficiency/");
	system("mv *.pdf ~/www/ttH/Full2016Dataset/Leptonic/Efficiency/");





	cout << endl << "Lufthansa, partner of Star Alliance, thanks you for choosing our company, we hope to see you again on board of our airctafts" << endl << endl;



























}









     






