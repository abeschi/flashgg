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


	allMVA -> SaveAs(title + ".pdf");
	allMVA -> SaveAs(title + ".png");




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

	RecoMVA -> SaveAs(title + "Reco.pdf");
	RecoMVA -> SaveAs(title + "Reco.png");



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

	TString SignalFold = "/afs/cern.ch/work/a/abeschi/ttHEfficiency_v9/";

	TChain* ttH;
	
	float lumiFactor = 35.9;

	ttH = new TChain("TTHSemiLeptonicEfficiencyDumper/trees/tth_13TeV_all");
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

	float fggmu_pt[2] = {0., 0.};
	float fggmu_eta[2] = {0., 0.};
	float fggmu_phi[2] = {0., 0.};
	float fggmu_isLoose[2] = {0., 0.};
	float fggmu_isMedium[2] = {0., 0.};
	float fggmu_MiniIso[2] = {0., 0.};
	float fggmu_trackIso[2] = {0., 0.};
	float fggmu_sumChargedHadronPt[2] = {0., 0.};
	float fggmu_sumNeutralHadronEt[2] = {0., 0.};
	float fggmu_sumPhotonEt[2] = {0., 0.};
	float fggmu_sumPUPt[2] = {0., 0.};
	float fggele_pt[2] = {0., 0.};
	float fggele_eta[2] = {0., 0.};
	float fggele_phi[2] = {0., 0.};
	float fggele_energy[2] = {0., 0.};
	float fggele_passLooseId[2] = {0., 0.};
	float fggele_passMediumId[2] = {0., 0.};
	float fggele_passTightId[2] = {0., 0.};
	float fggele_MVAMediumId[2] = {0., 0.};
	float fggele_MVATightId[2] = {0., 0.};
	float fggele_MiniIso[2] = {0., 0.};
	float fggele_ecalEnergy[2] = {0., 0.};
	float fggele_SCx[2] = {0., 0.};
	float fggele_SCy[2] = {0., 0.};
	float fggele_SCz[2] = {0., 0.};
	float fggele_SCeta[2] = {0., 0.};
	float fggele_SCphi[2] = {0., 0.};
	float fggele_dEtaSCTrackAtVtx[2] = {0., 0.};
	float fggele_dPhiSCTrackAtVtx[2] = {0., 0.};


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

		ttH -> SetBranchAddress(("fggmu_pt"+ std::to_string(i)).c_str(), &fggmu_pt[i-1]);
		ttH -> SetBranchAddress(("fggmu_eta"+ std::to_string(i)).c_str(), &fggmu_eta[i-1]);
		ttH -> SetBranchAddress(("fggmu_phi"+ std::to_string(i)).c_str(), &fggmu_phi[i-1]);
		ttH -> SetBranchAddress(("fggmu_isLoose"+ std::to_string(i)).c_str(), &fggmu_isLoose[i-1]);
		ttH -> SetBranchAddress(("fggmu_isMedium"+ std::to_string(i)).c_str(), &fggmu_isMedium[i-1]);
		ttH -> SetBranchAddress(("fggmu_MiniIso"+ std::to_string(i)).c_str(), &fggmu_MiniIso[i-1]);
		ttH -> SetBranchAddress(("fggmu_trackIso"+ std::to_string(i)).c_str(), &fggmu_trackIso[i-1]);
		ttH -> SetBranchAddress(("fggmu_sumChargedHadronPt"+ std::to_string(i)).c_str(), &fggmu_sumChargedHadronPt[i-1]);
		ttH -> SetBranchAddress(("fggmu_sumNeutralHadronEt"+ std::to_string(i)).c_str(), &fggmu_sumNeutralHadronEt[i-1]);
		ttH -> SetBranchAddress(("fggmu_sumPhotonEt"+ std::to_string(i)).c_str(), &fggmu_sumPhotonEt[i-1]);
		ttH -> SetBranchAddress(("fggmu_sumPUPt"+ std::to_string(i)).c_str(), &fggmu_sumPUPt[i-1]);

		ttH -> SetBranchAddress(("fggele_pt"+ std::to_string(i)).c_str(), &fggele_pt[i-1]);
		ttH -> SetBranchAddress(("fggele_eta"+ std::to_string(i)).c_str(), &fggele_eta[i-1]);
		ttH -> SetBranchAddress(("fggele_phi"+ std::to_string(i)).c_str(), &fggele_phi[i-1]);
		ttH -> SetBranchAddress(("fggele_energy"+ std::to_string(i)).c_str(), &fggele_energy[i-1]);
		ttH -> SetBranchAddress(("fggele_SCeta"+ std::to_string(i)).c_str(), &fggele_SCeta[i-1]);
		ttH -> SetBranchAddress(("fggele_SCphi"+ std::to_string(i)).c_str(), &fggele_SCphi[i-1]);
		ttH -> SetBranchAddress(("fggele_passLooseId"+ std::to_string(i)).c_str(), &fggele_passLooseId[i-1]);
		ttH -> SetBranchAddress(("fggele_passMediumId"+ std::to_string(i)).c_str(), &fggele_passMediumId[i-1]);
		ttH -> SetBranchAddress(("fggele_passTightId"+ std::to_string(i)).c_str(), &fggele_passTightId[i-1]);
		ttH -> SetBranchAddress(("fggele_MVAMediumId"+ std::to_string(i)).c_str(), &fggele_MVAMediumId[i-1]);
		ttH -> SetBranchAddress(("fggele_MVATightId"+ std::to_string(i)).c_str(), &fggele_MVATightId[i-1]);
		ttH -> SetBranchAddress(("fggele_MiniIso"+ std::to_string(i)).c_str(), &fggele_MiniIso[i-1]);
		ttH -> SetBranchAddress(("fggele_ecalEnergy"+ std::to_string(i)).c_str(), &fggele_ecalEnergy[i-1]);
		ttH -> SetBranchAddress(("fggele_SCx"+ std::to_string(i)).c_str(), &fggele_SCx[i-1]);
		ttH -> SetBranchAddress(("fggele_SCy"+ std::to_string(i)).c_str(), &fggele_SCy[i-1]);
		ttH -> SetBranchAddress(("fggele_SCz"+ std::to_string(i)).c_str(), &fggele_SCz[i-1]);
		ttH -> SetBranchAddress(("fggele_dEtaSCTrackAtVtx"+ std::to_string(i)).c_str(), &fggele_dEtaSCTrackAtVtx[i-1]);
		ttH -> SetBranchAddress(("fggele_dPhiSCTrackAtVtx"+ std::to_string(i)).c_str(), &fggele_dPhiSCTrackAtVtx[i-1]);
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

	const int n = 13;

	float efficiencyDiphotonVtx = 0.;
	float efficiencyVtx0 = 0.;
	float totalEvents = 0.;


//Study of muons

	float muTotal = 0;
	float muAcceptance = 0;
	float muMatched = 0;
	float muSelected = 0;
	float MuTest = 0;
	float muByStep[6] = {0., 0., 0., 0., 0., 0.};

	for(int ev=0; ev<nentries; ev++)
	{	
		ttH -> GetEntry(ev);
		if(ev%10000==0) cout << "Processing " << ev << "th event out of " << nentries << "\r" << flush;

		int selectedMuon = -1;

		if(abs(GenLepton_pdgId[0])!=13) continue;
		muTotal+=weight*lumiFactor;
		if(abs(GenLepton_eta[0])>2.4) continue;
		muAcceptance+=weight*lumiFactor;
		if(fggmu_pt[0]<0.) continue;
		float dr1 = compute_R(GenLepton_eta[0], fggmu_eta[0], GenLepton_phi[0], fggmu_phi[0]);
		float dr2 = compute_R(GenLepton_eta[0], fggmu_eta[1], GenLepton_phi[0], fggmu_phi[1]);
		float dr = min(dr1, dr2);
		if(dr==dr1) selectedMuon = 0;
		else selectedMuon = 1;
		if(dr>0.1) continue;
		muMatched+=weight*lumiFactor;

		float drReco = compute_R(lepton_eta[selectedMuon], fggmu_eta[selectedMuon], lepton_phi[selectedMuon], fggmu_phi[selectedMuon]);
		if(drReco<0.01) 
			muSelected+=weight*lumiFactor;

		if(fggmu_pt[selectedMuon]>20. && fggmu_MiniIso[selectedMuon]<0.06 && fggmu_isMedium[selectedMuon])
			MuTest+=weight*lumiFactor;

		if(fggmu_pt[selectedMuon]<20.) continue;
		muByStep[0]+=weight*lumiFactor;
		if(!fggmu_isLoose[selectedMuon]) continue;
		muByStep[1]+=weight*lumiFactor;
		if(!fggmu_isMedium[selectedMuon]) continue;
		muByStep[2]+=weight*lumiFactor;

		if(fggmu_MiniIso[selectedMuon]<0.06) 
			muByStep[3]+=weight*lumiFactor;
		if(fggmu_trackIso[selectedMuon]<0.10) 
			muByStep[4]+=weight*lumiFactor;

		float iso = fggmu_sumChargedHadronPt[selectedMuon] + max(0., fggmu_sumNeutralHadronEt[selectedMuon] + fggmu_sumPhotonEt[selectedMuon] - 0.5*fggmu_sumPUPt[selectedMuon]);
		if(iso/fggmu_pt[selectedMuon]<0.25) 
			muByStep[5]+=weight*lumiFactor;
	}

	cout << "Muons events: " << muTotal <<", within acceptance: " << muAcceptance << ", with matched reco candidate: " << muMatched << ", selected: " << muSelected << endl;
	cout << "Pt > 20:        " << muByStep[0] << endl;
	cout << "isLoose:        " << muByStep[1] << endl;
	cout << "isMedium:       " << muByStep[2] << endl;
	cout << "MiniIsolation:  " << muByStep[3] << endl;
	cout << "TrackIsolation: " << muByStep[4] << endl;
	cout << "StdIsolation:   " << muByStep[5] << endl;
	cout << MuTest << endl;



//Study of electrons

	float eleTotal = 0;
	float eleAcceptance = 0;
	float eleMatched = 0;
	float eleSelected = 0;
	float eleByStep[7] = {0., 0., 0., 0., 0., 0., 0.};

	for(int ev=0; ev<nentries; ev++)
	{	
		ttH -> GetEntry(ev);
		if(ev%10000==0) cout << "Processing " << ev << "th event out of " << nentries << "\r" << flush;

		int selectedEle = -1;

		if(abs(GenLepton_pdgId[0])!=11) continue;
		eleTotal+=weight*lumiFactor;
		if(abs(GenLepton_eta[0])>2.4 || (abs(GenLepton_eta[0])>1.4442 && abs(GenLepton_eta[0])<1.556) ) continue;
		eleAcceptance+=weight*lumiFactor;
		if(fggele_pt[0]<0.) continue;
		float dr1 = compute_R(GenLepton_eta[0], fggele_eta[0], GenLepton_phi[0], fggele_phi[0]);
		float dr2 = compute_R(GenLepton_eta[0], fggele_eta[1], GenLepton_phi[0], fggele_phi[1]);
		float dr = min(dr1, dr2);
		if(dr==dr1) selectedEle = 0;
		else selectedEle = 1;
		if(dr>0.1) continue;
		eleMatched+=weight*lumiFactor;

		float drReco = compute_R(lepton_eta[selectedEle+2], fggele_eta[selectedEle], lepton_phi[selectedEle+2], fggele_phi[selectedEle]);
		if(drReco<0.01) 
			eleSelected+=weight*lumiFactor;

		if(fggele_pt[selectedEle]<20.) continue;
		eleByStep[0]+=weight*lumiFactor;
		if(!fggele_passLooseId[selectedEle]) continue;
		eleByStep[1]+=weight*lumiFactor;
		if( (fabs(fggele_eta[selectedEle]) <= 1.479  && fggele_MiniIso[selectedEle] < 0.045) || (fabs(fggele_eta[selectedEle]) > 1.479  && fggele_MiniIso[selectedEle] < 0.08) );
			eleByStep[2]+=weight*lumiFactor;

		float drLead = compute_R(dipho_leadEta, fggele_eta[selectedEle], dipho_leadPhi, fggele_phi[selectedEle]);
		float drSublead = compute_R(dipho_subleadEta, fggele_eta[selectedEle], dipho_subleadPhi, fggele_phi[selectedEle]);
		if(drLead<0.35 || drSublead<0.35) continue;
		eleByStep[3]+=weight*lumiFactor;

		TLorentzVector ele;
		TLorentzVector eleSC;
		TLorentzVector leadPh;
		TLorentzVector subleadPh;

		ele.SetPtEtaPhiE(fggele_pt[selectedEle], fggele_eta[selectedEle], fggele_phi[selectedEle], fggele_energy[selectedEle]);
		eleSC.SetXYZT(fggele_SCx[selectedEle], fggele_SCy[selectedEle], fggele_SCz[selectedEle], fggele_ecalEnergy[selectedEle]);
		leadPh.SetPtEtaPhiE(dipho_leadPt, dipho_leadEta, dipho_leadPhi, dipho_leadEnergy);
		subleadPh.SetPtEtaPhiE(dipho_subleadPt, dipho_subleadEta, dipho_subleadPhi, dipho_subleadEnergy);

		if(leadPh.DeltaR(ele) < 0.35) continue;
		if(leadPh.DeltaR(eleSC)< 0.35) continue;
		if(subleadPh.DeltaR(ele)< 0.35) continue;
		if(subleadPh.DeltaR(eleSC)< 0.35) continue;
		eleByStep[4]+=weight*lumiFactor;

		float drTrack = sqrt(fggele_dEtaSCTrackAtVtx[selectedEle]*fggele_dEtaSCTrackAtVtx[selectedEle] + fggele_dPhiSCTrackAtVtx[selectedEle]*fggele_dPhiSCTrackAtVtx[selectedEle]);
		if(drTrack>0.35) continue;
		eleByStep[5]+=weight*lumiFactor;

		float M1 = (ele + leadPh).M();
		float M2 = (ele + leadPh).M();

		if(abs(M1-MZ)<5. || abs(M2-MZ)<5.) continue;
		eleByStep[6]+=weight*lumiFactor;
	}

	cout << "Electrons events: " << eleTotal <<", within acceptance: " << eleAcceptance << ", with matched reco candidate: " << eleMatched << ", selected: " << eleSelected << endl;
	cout << "Pt > 20:        " << eleByStep[0] << endl;
	cout << "isLoose:        " << eleByStep[1] << endl;
	cout << "MiniIsolation:  " << eleByStep[2] << endl;
	cout << "dr photons:     " << eleByStep[3] << endl;
	cout << "dr photons 2 :  " << eleByStep[4] << endl;
	cout << "dr tracks:      " << eleByStep[5] << endl;
	cout << "dM:             " << eleByStep[6] << endl;
//	cout << MuTest << endl;



//Study of taus

	float tauTotal = 0;
	float tauSelected = 0;

	for(int ev=0; ev<nentries; ev++)
	{	
		ttH -> GetEntry(ev);
		if(ev%10000==0) cout << "Processing " << ev << "th event out of " << nentries << "\r" << flush;

		int selectedEle = -1;

		if(abs(GenLepton_pdgId[0])!=15) continue;
		tauTotal+=weight*lumiFactor;
		if(lepton_pt[0]>0 ||  lepton_pt[2]>0)
			tauSelected+=weight*lumiFactor;

	}

	cout << "Tau events: " << eleTotal << ", selected: " << tauSelected << endl;

//	system("mv *.png ~/www/ttH/Full2016Dataset/Leptonic/Efficiency/");
//	system("mv *.pdf ~/www/ttH/Full2016Dataset/Leptonic/Efficiency/");





	cout << endl << "Lufthansa, partner of Star Alliance, thanks you for choosing our company, we hope to see you again on board of our airctafts" << endl << endl;



























}









     






